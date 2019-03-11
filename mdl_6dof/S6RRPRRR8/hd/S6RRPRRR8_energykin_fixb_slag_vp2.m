% Calculate kinetic energy for
% S6RRPRRR8
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:03
% EndTime: 2019-03-09 14:02:04
% DurationCPUTime: 0.66s
% Computational Cost: add. (1118->111), mult. (2493->172), div. (0->0), fcn. (1874->10), ass. (0->45)
t131 = pkin(7) * mrSges(3,3);
t120 = sin(qJ(5));
t124 = cos(qJ(5));
t117 = sin(pkin(11));
t118 = cos(pkin(11));
t122 = sin(qJ(2));
t130 = qJD(1) * t122;
t108 = qJD(2) * t118 - t117 * t130;
t109 = qJD(2) * t117 + t118 * t130;
t121 = sin(qJ(4));
t125 = cos(qJ(4));
t101 = t108 * t121 + t109 * t125;
t126 = cos(qJ(2));
t129 = qJD(1) * t126;
t115 = qJD(4) - t129;
t107 = (-pkin(2) * t126 - qJ(3) * t122 - pkin(1)) * qJD(1);
t113 = pkin(7) * t129 + qJD(2) * qJ(3);
t102 = t118 * t107 - t113 * t117;
t97 = -pkin(3) * t129 - pkin(8) * t109 + t102;
t103 = t117 * t107 + t118 * t113;
t99 = pkin(8) * t108 + t103;
t90 = -t121 * t99 + t125 * t97;
t86 = pkin(4) * t115 - pkin(9) * t101 + t90;
t100 = t108 * t125 - t109 * t121;
t91 = t121 * t97 + t125 * t99;
t88 = pkin(9) * t100 + t91;
t81 = t120 * t86 + t124 * t88;
t80 = -t120 * t88 + t124 * t86;
t111 = -qJD(2) * pkin(2) + pkin(7) * t130 + qJD(3);
t114 = qJD(5) + t115;
t104 = -pkin(3) * t108 + t111;
t94 = -pkin(4) * t100 + t104;
t123 = cos(qJ(6));
t119 = sin(qJ(6));
t112 = qJD(6) + t114;
t93 = t100 * t120 + t101 * t124;
t92 = t100 * t124 - t101 * t120;
t89 = -pkin(5) * t92 + t94;
t83 = t119 * t92 + t123 * t93;
t82 = -t119 * t93 + t123 * t92;
t79 = pkin(10) * t92 + t81;
t78 = pkin(5) * t114 - pkin(10) * t93 + t80;
t77 = t119 * t78 + t123 * t79;
t76 = -t119 * t79 + t123 * t78;
t1 = m(4) * (t102 ^ 2 + t103 ^ 2 + t111 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(5) * (t104 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t76 ^ 2 + t77 ^ 2 + t89 ^ 2) / 0.2e1 + (t94 * mrSges(6,2) - t80 * mrSges(6,3) + Ifges(6,1) * t93 / 0.2e1) * t93 + (t89 * mrSges(7,2) - t76 * mrSges(7,3) + Ifges(7,1) * t83 / 0.2e1) * t83 + (t90 * mrSges(5,1) - t91 * mrSges(5,2) + Ifges(5,3) * t115 / 0.2e1) * t115 + (t111 * mrSges(4,2) - t102 * mrSges(4,3) + Ifges(4,1) * t109 / 0.2e1) * t109 + (-t94 * mrSges(6,1) + t81 * mrSges(6,3) + Ifges(6,4) * t93 + Ifges(6,2) * t92 / 0.2e1) * t92 + (-t89 * mrSges(7,1) + t77 * mrSges(7,3) + Ifges(7,4) * t83 + Ifges(7,2) * t82 / 0.2e1) * t82 + (t104 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,5) * t115 + Ifges(5,1) * t101 / 0.2e1) * t101 + (t80 * mrSges(6,1) - t81 * mrSges(6,2) + Ifges(6,5) * t93 + Ifges(6,6) * t92 + Ifges(6,3) * t114 / 0.2e1) * t114 + (t76 * mrSges(7,1) - t77 * mrSges(7,2) + Ifges(7,5) * t83 + Ifges(7,6) * t82 + Ifges(7,3) * t112 / 0.2e1) * t112 + (-Ifges(4,6) * t129 - t111 * mrSges(4,1) + t103 * mrSges(4,3) + Ifges(4,4) * t109 + Ifges(4,2) * t108 / 0.2e1) * t108 + (-t104 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t101 + Ifges(5,6) * t115 + Ifges(5,2) * t100 / 0.2e1) * t100 + ((-pkin(7) * mrSges(3,1) + Ifges(3,5)) * qJD(2) * t122 + (-t102 * mrSges(4,1) + t103 * mrSges(4,2) - t109 * Ifges(4,5) + (-pkin(7) * mrSges(3,2) + Ifges(3,6)) * qJD(2)) * t126 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t122 ^ 2 + t126 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + t131) * t122) * t122 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t122 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + t131) * t126) * t126) * qJD(1)) * qJD(1);
T  = t1;
