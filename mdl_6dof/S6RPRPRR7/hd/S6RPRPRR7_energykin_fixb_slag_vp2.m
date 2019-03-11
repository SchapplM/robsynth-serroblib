% Calculate kinetic energy for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:48
% EndTime: 2019-03-09 03:54:48
% DurationCPUTime: 0.40s
% Computational Cost: add. (553->94), mult. (1134->148), div. (0->0), fcn. (752->8), ass. (0->40)
t97 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t114 = t97 * mrSges(4,3);
t113 = qJD(1) / 0.2e1;
t106 = sin(qJ(5));
t109 = cos(qJ(5));
t103 = sin(pkin(10));
t104 = cos(pkin(10));
t110 = cos(qJ(3));
t112 = -qJ(4) * qJD(1) + t97;
t91 = qJD(3) * pkin(3) + t112 * t110;
t107 = sin(qJ(3));
t92 = t112 * t107;
t82 = -t103 * t92 + t104 * t91;
t94 = (-t103 * t107 + t104 * t110) * qJD(1);
t78 = qJD(3) * pkin(4) - pkin(8) * t94 + t82;
t83 = t103 * t91 + t104 * t92;
t93 = (-t103 * t110 - t104 * t107) * qJD(1);
t79 = pkin(8) * t93 + t83;
t74 = t106 * t78 + t109 * t79;
t102 = qJD(1) * qJ(2);
t95 = t107 * qJD(1) * pkin(3) + qJD(4) + t102;
t73 = -t106 * t79 + t109 * t78;
t85 = -t106 * t94 + t109 * t93;
t87 = -pkin(4) * t93 + t95;
t111 = qJD(1) ^ 2;
t108 = cos(qJ(6));
t105 = sin(qJ(6));
t101 = t111 * qJ(2) ^ 2;
t100 = qJD(3) + qJD(5);
t99 = -qJD(1) * pkin(1) + qJD(2);
t86 = t106 * t93 + t109 * t94;
t84 = qJD(6) - t85;
t81 = t100 * t105 + t108 * t86;
t80 = t100 * t108 - t105 * t86;
t75 = -pkin(5) * t85 - pkin(9) * t86 + t87;
t72 = pkin(9) * t100 + t74;
t71 = -pkin(5) * t100 - t73;
t70 = t105 * t75 + t108 * t72;
t69 = -t105 * t72 + t108 * t75;
t1 = m(3) * (t99 ^ 2 + t101) / 0.2e1 + m(4) * (t101 + (t107 ^ 2 + t110 ^ 2) * t97 ^ 2) / 0.2e1 + m(5) * (t82 ^ 2 + t83 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t73 ^ 2 + t74 ^ 2 + t87 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + (t95 * mrSges(5,2) - t82 * mrSges(5,3) + Ifges(5,1) * t94 / 0.2e1) * t94 + (t87 * mrSges(6,2) - t73 * mrSges(6,3) + Ifges(6,1) * t86 / 0.2e1) * t86 + (t69 * mrSges(7,1) - t70 * mrSges(7,2) + Ifges(7,3) * t84 / 0.2e1) * t84 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t111 + (-t95 * mrSges(5,1) + t83 * mrSges(5,3) + Ifges(5,4) * t94 + Ifges(5,2) * t93 / 0.2e1) * t93 + (-t87 * mrSges(6,1) + t74 * mrSges(6,3) + Ifges(6,4) * t86 + Ifges(6,2) * t85 / 0.2e1) * t85 + (t71 * mrSges(7,2) - t69 * mrSges(7,3) + Ifges(7,5) * t84 + Ifges(7,1) * t81 / 0.2e1) * t81 + (-t71 * mrSges(7,1) + t70 * mrSges(7,3) + Ifges(7,4) * t81 + Ifges(7,6) * t84 + Ifges(7,2) * t80 / 0.2e1) * t80 + (t73 * mrSges(6,1) - t74 * mrSges(6,2) + Ifges(6,5) * t86 + Ifges(6,6) * t85 + Ifges(6,3) * t100 / 0.2e1) * t100 + (t99 * mrSges(3,2) + (mrSges(4,2) * t102 + (Ifges(4,1) * t113 - t114) * t110) * t110 + ((qJ(2) * mrSges(4,1) - Ifges(4,4) * t110) * qJD(1) + (Ifges(4,2) * t113 - t114) * t107) * t107) * qJD(1) + (t82 * mrSges(5,1) - t83 * mrSges(5,2) + Ifges(5,5) * t94 + Ifges(5,6) * t93 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (t110 * mrSges(4,1) - t107 * mrSges(4,2)) * t97 + (Ifges(4,5) * t110 - Ifges(4,6) * t107) * qJD(1)) * qJD(3);
T  = t1;
