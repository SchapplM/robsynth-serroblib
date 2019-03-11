% Calculate kinetic energy for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:21:40
% EndTime: 2019-03-09 13:21:41
% DurationCPUTime: 0.65s
% Computational Cost: add. (970->110), mult. (2231->171), div. (0->0), fcn. (1678->10), ass. (0->46)
t132 = pkin(7) * mrSges(3,3);
t131 = pkin(7) + qJ(3);
t120 = sin(qJ(5));
t124 = cos(qJ(5));
t117 = sin(pkin(11));
t118 = cos(pkin(11));
t122 = sin(qJ(2));
t126 = cos(qJ(2));
t110 = (t117 * t126 + t118 * t122) * qJD(1);
t121 = sin(qJ(4));
t125 = cos(qJ(4));
t104 = qJD(2) * t121 + t110 * t125;
t129 = qJD(1) * t126;
t130 = qJD(1) * t122;
t109 = -t117 * t130 + t118 * t129;
t108 = qJD(4) - t109;
t113 = qJD(2) * pkin(2) - t130 * t131;
t114 = t131 * t129;
t102 = t117 * t113 + t118 * t114;
t100 = qJD(2) * pkin(8) + t102;
t115 = qJD(3) + (-pkin(2) * t126 - pkin(1)) * qJD(1);
t97 = -pkin(3) * t109 - pkin(8) * t110 + t115;
t90 = -t100 * t121 + t125 * t97;
t84 = pkin(4) * t108 - pkin(9) * t104 + t90;
t103 = qJD(2) * t125 - t110 * t121;
t91 = t125 * t100 + t121 * t97;
t89 = pkin(9) * t103 + t91;
t81 = t120 * t84 + t124 * t89;
t80 = -t120 * t89 + t124 * t84;
t101 = t113 * t118 - t117 * t114;
t99 = -qJD(2) * pkin(3) - t101;
t106 = qJD(5) + t108;
t92 = -pkin(4) * t103 + t99;
t123 = cos(qJ(6));
t119 = sin(qJ(6));
t105 = qJD(6) + t106;
t94 = t103 * t120 + t104 * t124;
t93 = t103 * t124 - t104 * t120;
t88 = t119 * t93 + t123 * t94;
t87 = -t119 * t94 + t123 * t93;
t85 = -pkin(5) * t93 + t92;
t79 = pkin(10) * t93 + t81;
t78 = pkin(5) * t106 - pkin(10) * t94 + t80;
t77 = t119 * t78 + t123 * t79;
t76 = -t119 * t79 + t123 * t78;
t1 = m(4) * (t101 ^ 2 + t102 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t91 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t76 ^ 2 + t77 ^ 2 + t85 ^ 2) / 0.2e1 + (t92 * mrSges(6,2) - t80 * mrSges(6,3) + Ifges(6,1) * t94 / 0.2e1) * t94 + (t85 * mrSges(7,2) - t76 * mrSges(7,3) + Ifges(7,1) * t88 / 0.2e1) * t88 + (t115 * mrSges(4,2) - t101 * mrSges(4,3) + Ifges(4,1) * t110 / 0.2e1) * t110 + (t90 * mrSges(5,1) - t91 * mrSges(5,2) + Ifges(5,3) * t108 / 0.2e1) * t108 + (-t92 * mrSges(6,1) + t81 * mrSges(6,3) + Ifges(6,4) * t94 + Ifges(6,2) * t93 / 0.2e1) * t93 + (-t85 * mrSges(7,1) + t77 * mrSges(7,3) + Ifges(7,4) * t88 + Ifges(7,2) * t87 / 0.2e1) * t87 + (-t115 * mrSges(4,1) + t102 * mrSges(4,3) + Ifges(4,4) * t110 + Ifges(4,2) * t109 / 0.2e1) * t109 + (t99 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,5) * t108 + Ifges(5,1) * t104 / 0.2e1) * t104 + (t80 * mrSges(6,1) - t81 * mrSges(6,2) + Ifges(6,5) * t94 + Ifges(6,6) * t93 + Ifges(6,3) * t106 / 0.2e1) * t106 + (t76 * mrSges(7,1) - t77 * mrSges(7,2) + Ifges(7,5) * t88 + Ifges(7,6) * t87 + Ifges(7,3) * t105 / 0.2e1) * t105 + (-t99 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t104 + Ifges(5,6) * t108 + Ifges(5,2) * t103 / 0.2e1) * t103 + (t101 * mrSges(4,1) - t102 * mrSges(4,2) + Ifges(4,5) * t110 + Ifges(4,6) * t109 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t122 + Ifges(3,6) * t126 + (-mrSges(3,1) * t122 - mrSges(3,2) * t126) * pkin(7)) * qJD(1)) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t122 ^ 2 + t126 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t132) * t126) * t126 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t126 + (Ifges(3,1) / 0.2e1 + t132) * t122) * t122) * qJD(1) ^ 2;
T  = t1;
