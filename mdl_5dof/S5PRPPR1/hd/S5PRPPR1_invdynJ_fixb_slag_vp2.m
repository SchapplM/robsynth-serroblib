% Calculate vector of inverse dynamics joint torques for
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:47
% EndTime: 2019-12-05 15:21:56
% DurationCPUTime: 3.44s
% Computational Cost: add. (1462->266), mult. (3335->382), div. (0->0), fcn. (2283->10), ass. (0->130)
t104 = sin(pkin(8));
t106 = cos(pkin(8));
t175 = t104 ^ 2 + t106 ^ 2;
t85 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t157 = m(5) + m(6);
t137 = m(4) + t157;
t102 = pkin(7) + qJ(2);
t98 = cos(t102);
t174 = g(2) * t98;
t103 = sin(pkin(9));
t108 = sin(qJ(5));
t105 = cos(pkin(9));
t109 = cos(qJ(5));
t143 = t105 * t109;
t118 = t103 * t108 - t143;
t74 = t103 * t109 + t105 * t108;
t66 = t74 * qJD(5);
t23 = (-qJD(2) * t66 - qJDD(2) * t118) * t104;
t173 = Ifges(6,5) * t23;
t166 = t118 * qJD(5);
t24 = (qJD(2) * t166 - qJDD(2) * t74) * t104;
t172 = Ifges(6,6) * t24;
t133 = qJDD(2) * t106;
t86 = qJDD(5) - t133;
t171 = Ifges(6,3) * t86;
t141 = qJD(2) * t106;
t170 = t118 * t141 - t166;
t113 = qJD(2) * t74;
t169 = t106 * t113 - t66;
t59 = qJDD(1) * t106 - t104 * t85;
t60 = qJDD(1) * t104 + t106 * t85;
t167 = -t104 * t59 + t106 * t60;
t122 = mrSges(4,1) * t106 - mrSges(4,2) * t104;
t164 = t104 * mrSges(6,3) + mrSges(3,1) + t122;
t121 = t103 * mrSges(5,1) + t105 * mrSges(5,2);
t114 = t121 * t104;
t142 = qJD(2) * t104;
t128 = t103 * t142;
t76 = -qJ(3) * t142 + qJD(1) * t106;
t72 = qJD(4) - t76;
t43 = pkin(4) * t128 + t72;
t44 = t104 * t113;
t46 = -t108 * t128 + t142 * t143;
t163 = m(6) * t43 + mrSges(6,1) * t44 + mrSges(6,2) * t46 + qJD(2) * t114;
t145 = t104 * t105;
t131 = pkin(6) * t145;
t117 = -pkin(4) * t106 - t131;
t138 = qJD(4) * t104;
t125 = -qJ(4) * t104 - pkin(2);
t79 = -pkin(3) * t106 + t125;
t39 = -qJD(2) * t138 + qJDD(2) * t79 + qJDD(3);
t17 = -t103 * t60 + t105 * t39;
t10 = qJDD(2) * t117 + t17;
t134 = qJDD(2) * t104;
t127 = t103 * t134;
t18 = t103 * t39 + t105 * t60;
t11 = -pkin(6) * t127 + t18;
t58 = qJD(2) * t79 + qJD(3);
t77 = qJ(3) * t141 + qJD(1) * t104;
t26 = -t103 * t77 + t105 * t58;
t16 = qJD(2) * t117 + t26;
t27 = t103 * t58 + t105 * t77;
t25 = -pkin(6) * t128 + t27;
t5 = -t108 * t25 + t109 * t16;
t1 = qJD(5) * t5 + t10 * t108 + t109 * t11;
t6 = t108 * t16 + t109 * t25;
t2 = -qJD(5) * t6 + t10 * t109 - t108 * t11;
t162 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t154 = pkin(4) * t103;
t161 = -m(6) * t154 + mrSges(3,2) - mrSges(4,3) - t121;
t159 = t46 / 0.2e1;
t158 = m(2) + m(3);
t126 = t105 * t134;
t53 = mrSges(5,1) * t127 + mrSges(5,2) * t126;
t7 = -t24 * mrSges(6,1) + t23 * mrSges(6,2);
t156 = -t53 - t7;
t155 = Ifges(6,4) * t46;
t144 = t105 * t106;
t42 = qJ(3) * t144 + t103 * t79;
t96 = sin(t102);
t149 = t106 * t96;
t148 = t106 * t98;
t147 = t103 * t104;
t146 = t103 * t106;
t139 = qJD(3) * t106;
t132 = t171 + t172 + t173;
t123 = -mrSges(4,1) * t133 + mrSges(4,2) * t134;
t120 = -t104 * t76 + t106 * t77;
t71 = t105 * t79;
t28 = -t131 + t71 + (-qJ(3) * t103 - pkin(4)) * t106;
t31 = -pkin(6) * t147 + t42;
t8 = -t108 * t31 + t109 * t28;
t9 = t108 * t28 + t109 * t31;
t119 = -t104 * (-pkin(6) - qJ(4)) + t106 * (pkin(4) * t105 + pkin(3));
t56 = qJDD(4) - t59;
t116 = -mrSges(5,1) * t106 - mrSges(5,3) * t145;
t115 = mrSges(5,2) * t106 - mrSges(5,3) * t147;
t111 = m(5) * pkin(3) + t105 * mrSges(5,1) - t103 * mrSges(5,2);
t101 = pkin(9) + qJ(5);
t97 = cos(t101);
t95 = sin(t101);
t94 = -qJDD(2) * pkin(2) + qJDD(3);
t87 = qJD(5) - t141;
t75 = (qJ(3) + t154) * t104;
t68 = t116 * qJD(2);
t67 = t115 * qJD(2);
t64 = t116 * qJDD(2);
t63 = t115 * qJDD(2);
t62 = -t103 * t138 + t105 * t139;
t61 = -t103 * t139 - t105 * t138;
t55 = t118 * t104;
t54 = t74 * t104;
t49 = t104 * t166;
t48 = t104 * t66;
t41 = -qJ(3) * t146 + t71;
t40 = Ifges(6,4) * t44;
t38 = t148 * t97 + t95 * t96;
t37 = -t148 * t95 + t96 * t97;
t36 = -t149 * t97 + t95 * t98;
t35 = t149 * t95 + t97 * t98;
t34 = pkin(4) * t127 + t56;
t30 = mrSges(6,1) * t87 - mrSges(6,3) * t46;
t29 = -mrSges(6,2) * t87 - mrSges(6,3) * t44;
t15 = Ifges(6,1) * t46 + Ifges(6,5) * t87 - t40;
t14 = -Ifges(6,2) * t44 + Ifges(6,6) * t87 + t155;
t13 = -mrSges(6,2) * t86 + mrSges(6,3) * t24;
t12 = mrSges(6,1) * t86 - mrSges(6,3) * t23;
t4 = -qJD(5) * t9 - t108 * t62 + t109 * t61;
t3 = qJD(5) * t8 + t108 * t61 + t109 * t62;
t19 = [-t54 * t12 - t55 * t13 - t48 * t29 + t49 * t30 + t156 * t106 + (-t103 * t64 + t105 * t63) * t104 + t158 * qJDD(1) + m(4) * (t104 * t60 + t106 * t59) + m(5) * (-t106 * t56 + (-t103 * t17 + t105 * t18) * t104) + m(6) * (-t1 * t55 - t106 * t34 - t2 * t54 - t48 * t6 + t49 * t5) + (-t137 - t158) * g(3); m(4) * (-pkin(2) * t94 + qJ(3) * t167 + t120 * qJD(3)) + t163 * qJD(3) * t104 - (mrSges(6,2) * t34 - mrSges(6,3) * t2 + Ifges(6,1) * t23 + Ifges(6,4) * t24 + Ifges(6,5) * t86) * t55 + m(5) * (t17 * t41 + t18 * t42 + t26 * t61 + t27 * t62) + m(6) * (t1 * t9 + t2 * t8 + t3 * t6 + t34 * t75 + t4 * t5) + (-t38 * mrSges(6,1) - t37 * mrSges(6,2) + (-m(6) * t119 - t164) * t98 - t137 * (t98 * pkin(2) + t96 * qJ(3)) + t161 * t96) * g(2) + (-t36 * mrSges(6,1) - t35 * mrSges(6,2) + (-m(5) * t125 + t104 * mrSges(5,3) + t111 * t106 + m(4) * pkin(2) - m(6) * (-pkin(2) - t119) + t164) * t96 + (-qJ(3) * t137 + t161) * t98) * g(1) + (t175 * t85 + t167) * mrSges(4,3) + t56 * t114 + (t5 * t48 + t6 * t49) * mrSges(6,3) + (-Ifges(6,1) * t48 + Ifges(6,4) * t49) * t159 + t87 * (-Ifges(6,5) * t48 + Ifges(6,6) * t49) / 0.2e1 + t43 * (-mrSges(6,1) * t49 - mrSges(6,2) * t48) - t44 * (-Ifges(6,4) * t48 + Ifges(6,2) * t49) / 0.2e1 + t42 * t63 + t41 * t64 + t62 * t67 + t61 * t68 + t75 * t7 + t49 * t14 / 0.2e1 - pkin(2) * t123 - t94 * t122 + t18 * t115 + t17 * t116 - (-mrSges(6,1) * t34 + mrSges(6,3) * t1 + Ifges(6,4) * t23 + Ifges(6,2) * t24 + Ifges(6,6) * t86) * t54 + ((Ifges(5,1) * t105 - Ifges(5,4) * t103) * t126 + m(5) * (qJ(3) * t56 + qJD(3) * t72) + Ifges(4,1) * t134 - (Ifges(5,4) * t105 - Ifges(5,2) * t103) * t127 + qJ(3) * t53 - (m(5) * qJ(4) + mrSges(5,3)) * t174 + (-Ifges(5,5) * t105 + Ifges(5,6) * t103 + Ifges(4,4)) * t133) * t104 + (-Ifges(5,5) * t126 + Ifges(4,4) * t134 - t171 / 0.2e1 - t172 / 0.2e1 + Ifges(5,6) * t127 - t173 / 0.2e1 - t111 * t174 - t132 / 0.2e1 + (Ifges(4,2) + Ifges(5,3)) * t133 + t162) * t106 - t48 * t15 / 0.2e1 + t3 * t29 + t4 * t30 + t8 * t12 + t9 * t13 + Ifges(3,3) * qJDD(2); t103 * t63 + t105 * t64 - t118 * t12 + t74 * t13 + t169 * t30 + t170 * t29 + m(5) * (t103 * t18 + t105 * t17) + m(4) * t94 + t123 + (-g(1) * t96 + t174) * t137 + (t1 * t74 - t118 * t2 + t169 * t5 + t170 * t6) * m(6) + ((t103 * t68 - t105 * t67) * t106 - m(5) * (t144 * t27 - t146 * t26) - m(4) * t120 - t175 * mrSges(4,3) * qJD(2) + (-m(5) * t72 - t163) * t104) * qJD(2); t44 * t29 + t46 * t30 + t157 * t106 * g(3) + m(5) * t56 + ((t103 * t67 + t105 * t68 - m(5) * (-t103 * t27 - t105 * t26)) * qJD(2) + t157 * (-g(1) * t98 - g(2) * t96)) * t104 - t156 + (t44 * t6 + t46 * t5 + t34) * m(6); -t43 * (mrSges(6,1) * t46 - mrSges(6,2) * t44) - t46 * (-Ifges(6,1) * t44 - t155) / 0.2e1 + t14 * t159 - t87 * (-Ifges(6,5) * t44 - Ifges(6,6) * t46) / 0.2e1 - t5 * t29 + t6 * t30 - g(1) * (mrSges(6,1) * t37 - mrSges(6,2) * t38) - g(2) * (-mrSges(6,1) * t35 + mrSges(6,2) * t36) - g(3) * (-mrSges(6,1) * t95 - mrSges(6,2) * t97) * t104 + (-t44 * t5 + t46 * t6) * mrSges(6,3) + t132 + (-Ifges(6,2) * t46 + t15 - t40) * t44 / 0.2e1 - t162;];
tau = t19;
