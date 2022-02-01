% Calculate vector of inverse dynamics joint torques for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:28
% EndTime: 2022-01-23 09:34:32
% DurationCPUTime: 1.55s
% Computational Cost: add. (2691->241), mult. (4880->315), div. (0->0), fcn. (2672->16), ass. (0->118)
t103 = sin(qJ(5));
t159 = t103 / 0.2e1;
t170 = mrSges(5,2) - mrSges(6,3);
t107 = cos(qJ(5));
t138 = qJD(5) * t103;
t98 = qJDD(1) + qJDD(3);
t92 = qJDD(4) + t98;
t99 = qJD(1) + qJD(3);
t95 = qJD(4) + t99;
t55 = t107 * t92 - t138 * t95;
t137 = qJD(5) * t107;
t56 = t103 * t92 + t137 * t95;
t33 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t104 = sin(qJ(4));
t108 = cos(qJ(4));
t105 = sin(qJ(3));
t109 = cos(qJ(3));
t101 = sin(pkin(9));
t157 = pkin(1) * t101;
t132 = qJD(1) * t157;
t102 = cos(pkin(9));
t87 = pkin(1) * t102 + pkin(2);
t69 = t87 * qJD(1);
t47 = t105 * t69 + t109 * t132;
t141 = t108 * t47;
t46 = -t105 * t132 + t109 * t69;
t39 = pkin(3) * t99 + t46;
t18 = t104 * t39 + t141;
t135 = qJD(1) * qJD(3);
t140 = qJD(3) * t69;
t68 = t87 * qJDD(1);
t28 = -t105 * t140 + t109 * t68 + (-qJDD(1) * t105 - t109 * t135) * t157;
t20 = pkin(3) * t98 + t28;
t27 = t109 * t140 + t105 * t68 + (qJDD(1) * t109 - t105 * t135) * t157;
t9 = -qJD(4) * t18 - t104 * t27 + t108 * t20;
t6 = -pkin(4) * t92 - t9;
t169 = m(6) * t6 + t33;
t146 = t104 * t47;
t17 = t108 * t39 - t146;
t70 = -t107 * mrSges(6,1) + t103 * mrSges(6,2);
t168 = -mrSges(5,1) + t70;
t127 = mrSges(6,1) * t103 + mrSges(6,2) * t107;
t15 = -pkin(4) * t95 - t17;
t167 = t15 * t127 + qJD(5) * (Ifges(6,5) * t107 - Ifges(6,6) * t103) / 0.2e1;
t166 = m(3) * pkin(1);
t165 = t28 * mrSges(4,1) + Ifges(4,3) * t98;
t16 = pkin(8) * t95 + t18;
t14 = qJD(2) * t103 + t107 * t16;
t139 = t14 * qJD(5);
t8 = t17 * qJD(4) + t104 * t20 + t108 * t27;
t5 = pkin(8) * t92 + t8;
t3 = qJDD(2) * t107 - t103 * t5 - t139;
t164 = -t3 - t139;
t58 = -t105 * t157 + t109 * t87;
t13 = qJD(2) * t107 - t103 * t16;
t2 = qJD(5) * t13 + qJDD(2) * t103 + t107 * t5;
t156 = t107 * t2;
t163 = -t13 * t137 + t156;
t144 = t107 * t14;
t123 = -t103 * t13 + t144;
t161 = m(5) + m(6);
t100 = qJ(1) + pkin(9);
t96 = qJ(3) + t100;
t85 = sin(t96);
t160 = g(1) * t85;
t155 = t3 * t103;
t57 = pkin(3) + t58;
t59 = t105 * t87 + t109 * t157;
t32 = t104 * t57 + t108 * t59;
t88 = qJ(4) + t96;
t81 = sin(t88);
t82 = cos(t88);
t154 = t82 * pkin(4) + t81 * pkin(8);
t110 = cos(qJ(1));
t94 = cos(t100);
t153 = t110 * pkin(1) + pkin(2) * t94;
t152 = Ifges(6,4) * t103;
t151 = Ifges(6,4) * t107;
t150 = Ifges(6,2) * t107;
t147 = t103 * t95;
t143 = t107 * t95;
t136 = m(4) + t161;
t86 = cos(t96);
t80 = pkin(3) * t86;
t134 = t80 + t153;
t131 = m(3) + t136;
t54 = t70 * t95;
t129 = t95 * mrSges(5,1) - t54;
t126 = t150 + t152;
t124 = t103 * t14 + t107 * t13;
t31 = -t104 * t59 + t108 * t57;
t120 = t103 * (Ifges(6,1) * t107 - t152);
t62 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t147;
t63 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t143;
t119 = t95 * mrSges(5,2) + t103 * t62 - t107 * t63;
t118 = t168 * t82 + t170 * t81;
t117 = -t86 * mrSges(4,1) + t85 * mrSges(4,2) + t118;
t44 = Ifges(6,6) * qJD(5) + t126 * t95;
t72 = Ifges(6,4) * t143;
t45 = Ifges(6,1) * t147 + Ifges(6,5) * qJD(5) + t72;
t116 = t9 * mrSges(5,1) + mrSges(6,3) * t156 + t6 * t70 + (Ifges(6,1) * t56 + Ifges(6,4) * t55) * t159 + t107 * (Ifges(6,4) * t56 + Ifges(6,2) * t55) / 0.2e1 + t55 * t126 / 0.2e1 + t56 * (Ifges(6,1) * t103 + t151) / 0.2e1 - t44 * t138 / 0.2e1 + Ifges(5,3) * t92 + (t45 + t95 * (-Ifges(6,2) * t103 + t151)) * t137 / 0.2e1 + (0.2e1 * Ifges(6,5) * t159 + Ifges(6,6) * t107) * qJDD(5) + (t120 * t95 / 0.2e1 + t167) * qJD(5);
t115 = (m(6) * pkin(4) - t168) * t81 + (-m(6) * pkin(8) + t170) * t82;
t114 = t86 * mrSges(4,2) + t115;
t41 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t55;
t42 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t56;
t113 = -t62 * t137 - t63 * t138 + m(6) * (-t138 * t14 - t155 + t163) + t107 * t41 - t103 * t42;
t112 = (-qJD(5) * t124 - t155) * mrSges(6,3) + t116 - t8 * mrSges(5,2);
t106 = sin(qJ(1));
t93 = sin(t100);
t53 = t59 * qJD(3);
t52 = t58 * qJD(3);
t30 = pkin(8) + t32;
t29 = -pkin(4) - t31;
t22 = t108 * t46 - t146;
t21 = t104 * t46 + t141;
t11 = qJD(4) * t32 + t104 * t52 + t108 * t53;
t10 = qJD(4) * t31 - t104 * t53 + t108 * t52;
t1 = [(mrSges(2,2) * t106 - mrSges(3,1) * t94 + mrSges(3,2) * t93 - m(6) * (t134 + t154) - m(5) * t134 - m(4) * t153 + (-mrSges(2,1) - t166) * t110 + t117) * g(2) + t116 + (t10 * t63 + t30 * t41 + (-t13 * mrSges(6,3) - t30 * t62) * qJD(5)) * t107 + m(6) * (t10 * t144 + t11 * t15 + t163 * t30 + t29 * t6) + t11 * t54 + t29 * t33 + m(4) * (t27 * t59 + t28 * t58 - t46 * t53 + t47 * t52) + m(5) * (t10 * t18 - t11 * t17 + t31 * t9 + t32 * t8) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t102 * mrSges(3,1) - 0.2e1 * t101 * mrSges(3,2) + (t101 ^ 2 + t102 ^ 2) * t166) * pkin(1)) * qJDD(1) + (t164 * mrSges(6,3) + (-m(6) * t13 - t62) * t10 + (m(6) * t164 - qJD(5) * t63 - t42) * t30) * t103 + (-t11 * t95 + t31 * t92) * mrSges(5,1) + (mrSges(2,2) * t110 + mrSges(3,2) * t94 + (pkin(3) * t161 + mrSges(4,1)) * t85 + (pkin(2) * t136 + mrSges(3,1)) * t93 + (pkin(1) * t131 + mrSges(2,1)) * t106 + t114) * g(1) + (-t10 * t95 - t32 * t92 - t8) * mrSges(5,2) + (-t52 * t99 - t59 * t98 - t27) * mrSges(4,2) + (-t53 * t99 + t58 * t98) * mrSges(4,1) + t165; t63 * t137 - t62 * t138 + t103 * t41 + t107 * t42 + m(6) * (qJD(5) * t123 + t103 * t2 + t107 * t3) + (m(3) + m(4) + m(5)) * qJDD(2) - t131 * g(3); -m(6) * (t123 * t22 + t15 * t21) + t112 + (-m(6) * (t80 + t154) + t117) * g(2) + t47 * t99 * mrSges(4,1) - m(5) * (-t17 * t21 + t18 * t22) + t113 * (pkin(3) * t104 + pkin(8)) + (m(6) * t160 + (mrSges(5,1) * t108 - mrSges(5,2) * t104) * t92 + (-g(2) * t86 + t104 * t8 + t108 * t9 + t160) * m(5) + ((-m(5) * t17 + m(6) * t15 - t129) * t104 + (m(5) * t18 + m(6) * t123 - t119) * t108) * qJD(4)) * pkin(3) + (t85 * mrSges(4,1) + t114) * g(1) + (t46 * t99 - t27) * mrSges(4,2) + t119 * t22 + t129 * t21 + t165 + t169 * (-pkin(3) * t108 - pkin(4)); -m(6) * (t123 * t17 + t15 * t18) + t115 * g(1) + t129 * t18 + t119 * t17 + t112 + (-m(6) * t154 + t118) * g(2) + t113 * pkin(8) - t169 * pkin(4); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t56 + Ifges(6,6) * t55 + Ifges(6,3) * qJDD(5) + g(3) * t70 - t13 * t63 + t14 * t62 + (t44 * t159 + (-t120 / 0.2e1 + t150 * t159) * t95 + t124 * mrSges(6,3) - (t45 + t72) * t107 / 0.2e1 - t167) * t95 + (g(1) * t82 + g(2) * t81) * t127;];
tau = t1;
