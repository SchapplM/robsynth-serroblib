% Calculate kinetic energy for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:19
% EndTime: 2019-03-09 15:28:19
% DurationCPUTime: 0.46s
% Computational Cost: add. (424->105), mult. (877->138), div. (0->0), fcn. (542->6), ass. (0->37)
t114 = -pkin(4) - pkin(9);
t113 = -pkin(8) - pkin(7);
t112 = pkin(7) * mrSges(3,3);
t100 = sin(qJ(3));
t103 = cos(qJ(3));
t101 = sin(qJ(2));
t111 = qJD(1) * t101;
t91 = qJD(2) * pkin(2) + t113 * t111;
t104 = cos(qJ(2));
t110 = qJD(1) * t104;
t92 = t113 * t110;
t80 = t100 * t91 - t103 * t92;
t93 = -qJD(1) * pkin(1) - pkin(2) * t110;
t97 = qJD(2) + qJD(3);
t78 = t97 * qJ(4) + t80;
t79 = t100 * t92 + t103 * t91;
t87 = t100 * t111 - t103 * t110;
t88 = (t100 * t104 + t101 * t103) * qJD(1);
t76 = t87 * pkin(3) - t88 * qJ(4) + t93;
t109 = qJD(4) - t79;
t108 = qJD(5) - t76;
t75 = -qJ(5) * t87 - t78;
t107 = -qJ(5) * t88 + t109;
t102 = cos(qJ(6));
t99 = sin(qJ(6));
t86 = qJD(6) + t88;
t82 = t102 * t87 - t97 * t99;
t81 = -t102 * t97 - t87 * t99;
t77 = -pkin(3) * t97 + t109;
t74 = pkin(5) * t97 - t75;
t73 = -pkin(4) * t87 + t108;
t72 = (-pkin(3) - pkin(4)) * t97 + t107;
t71 = (-pkin(3) + t114) * t97 + t107;
t70 = pkin(5) * t88 + t114 * t87 + t108;
t69 = t102 * t71 + t70 * t99;
t68 = t102 * t70 - t71 * t99;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t79 ^ 2 + t80 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t68 ^ 2 + t69 ^ 2 + t74 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t75 ^ 2) / 0.2e1 + m(5) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + (t68 * mrSges(7,1) - t69 * mrSges(7,2) + Ifges(7,3) * t86 / 0.2e1) * t86 + (t74 * mrSges(7,2) - t68 * mrSges(7,3) + Ifges(7,5) * t86 + Ifges(7,1) * t82 / 0.2e1) * t82 + (-t74 * mrSges(7,1) + t69 * mrSges(7,3) + Ifges(7,4) * t82 + Ifges(7,6) * t86 + Ifges(7,2) * t81 / 0.2e1) * t81 + (t79 * mrSges(4,1) - t77 * mrSges(5,1) - t75 * mrSges(6,1) - t80 * mrSges(4,2) + t72 * mrSges(6,2) + t78 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t97) * t97 + (t73 * mrSges(6,1) + t93 * mrSges(4,2) + t77 * mrSges(5,2) - t79 * mrSges(4,3) - t76 * mrSges(5,3) - t72 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t88 + (Ifges(5,4) + Ifges(4,5) + Ifges(6,6)) * t97) * t88 + (t93 * mrSges(4,1) + t76 * mrSges(5,1) - t78 * mrSges(5,2) + t73 * mrSges(6,2) - t80 * mrSges(4,3) - t75 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t87 + (-Ifges(6,5) - Ifges(4,6) + Ifges(5,6)) * t97 + (-Ifges(4,4) - Ifges(6,4) + Ifges(5,5)) * t88) * t87 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t101 ^ 2 + t104 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t112 + Ifges(3,2) / 0.2e1) * t104) * t104 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t104 + (t112 + Ifges(3,1) / 0.2e1) * t101) * t101) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t104 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t101) * qJD(2)) * qJD(1);
T  = t1;
