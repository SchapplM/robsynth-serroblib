% Calculate kinetic energy for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:43
% EndTime: 2019-03-09 08:45:43
% DurationCPUTime: 0.54s
% Computational Cost: add. (560->107), mult. (1299->154), div. (0->0), fcn. (890->8), ass. (0->40)
t124 = pkin(7) * mrSges(3,3);
t123 = pkin(7) + qJ(3);
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t113 = sin(qJ(2));
t121 = qJD(1) * t113;
t102 = qJD(2) * pkin(2) - t121 * t123;
t116 = cos(qJ(2));
t120 = qJD(1) * t116;
t103 = t123 * t120;
t109 = sin(pkin(10));
t122 = cos(pkin(10));
t93 = t102 * t122 - t109 * t103;
t119 = qJD(4) - t93;
t99 = (t109 * t116 + t113 * t122) * qJD(1);
t82 = -t99 * pkin(8) + (-pkin(3) - pkin(4)) * qJD(2) + t119;
t94 = t109 * t102 + t122 * t103;
t92 = qJD(2) * qJ(4) + t94;
t98 = t109 * t121 - t120 * t122;
t84 = pkin(8) * t98 + t92;
t79 = t112 * t82 + t115 * t84;
t104 = -qJD(1) * pkin(1) - pkin(2) * t120 + qJD(3);
t87 = t98 * pkin(3) - t99 * qJ(4) + t104;
t78 = -t112 * t84 + t115 * t82;
t90 = -t112 * t99 + t115 * t98;
t80 = -pkin(4) * t98 - t87;
t114 = cos(qJ(6));
t111 = sin(qJ(6));
t107 = -qJD(2) + qJD(5);
t91 = t112 * t98 + t115 * t99;
t89 = -qJD(2) * pkin(3) + t119;
t88 = qJD(6) - t90;
t86 = t107 * t111 + t114 * t91;
t85 = t107 * t114 - t111 * t91;
t77 = pkin(9) * t107 + t79;
t76 = -pkin(5) * t107 - t78;
t75 = -pkin(5) * t90 - pkin(9) * t91 + t80;
t74 = t111 * t75 + t114 * t77;
t73 = -t111 * t77 + t114 * t75;
t1 = m(4) * (t104 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t87 ^ 2 + t89 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(7) * (t73 ^ 2 + t74 ^ 2 + t76 ^ 2) / 0.2e1 + (t80 * mrSges(6,2) - t78 * mrSges(6,3) + Ifges(6,1) * t91 / 0.2e1) * t91 + (t73 * mrSges(7,1) - t74 * mrSges(7,2) + Ifges(7,3) * t88 / 0.2e1) * t88 + (-t80 * mrSges(6,1) + t79 * mrSges(6,3) + Ifges(6,4) * t91 + Ifges(6,2) * t90 / 0.2e1) * t90 + (t76 * mrSges(7,2) - t73 * mrSges(7,3) + Ifges(7,5) * t88 + Ifges(7,1) * t86 / 0.2e1) * t86 + (-t76 * mrSges(7,1) + t74 * mrSges(7,3) + Ifges(7,4) * t86 + Ifges(7,6) * t88 + Ifges(7,2) * t85 / 0.2e1) * t85 + (t78 * mrSges(6,1) - t79 * mrSges(6,2) + Ifges(6,5) * t91 + Ifges(6,6) * t90 + Ifges(6,3) * t107 / 0.2e1) * t107 + (t104 * mrSges(4,2) + t89 * mrSges(5,2) - t93 * mrSges(4,3) - t87 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t99) * t99 + (t104 * mrSges(4,1) + t87 * mrSges(5,1) - t92 * mrSges(5,2) - t94 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t98 + (-Ifges(4,4) + Ifges(5,5)) * t99) * t98 + (t93 * mrSges(4,1) - t89 * mrSges(5,1) - t94 * mrSges(4,2) + t92 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * qJD(2) + (Ifges(5,4) + Ifges(4,5)) * t99 + (-Ifges(4,6) + Ifges(5,6)) * t98 + (Ifges(3,5) * t113 + Ifges(3,6) * t116 + (-mrSges(3,1) * t113 - mrSges(3,2) * t116) * pkin(7)) * qJD(1)) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t113 ^ 2 + t116 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t124 + Ifges(3,2) / 0.2e1) * t116) * t116 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t116 + (t124 + Ifges(3,1) / 0.2e1) * t113) * t113) * qJD(1) ^ 2;
T  = t1;
