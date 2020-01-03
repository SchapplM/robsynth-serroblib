% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:57
% EndTime: 2019-12-31 18:32:04
% DurationCPUTime: 3.00s
% Computational Cost: add. (2569->344), mult. (6821->462), div. (0->0), fcn. (4676->6), ass. (0->158)
t158 = mrSges(5,2) - mrSges(4,1);
t189 = mrSges(4,3) + mrSges(5,1);
t103 = sin(qJ(3));
t101 = cos(pkin(8));
t156 = pkin(6) + qJ(2);
t91 = t156 * t101;
t88 = qJD(1) * t91;
t141 = t103 * t88;
t167 = cos(qJ(3));
t100 = sin(pkin(8));
t90 = t156 * t100;
t87 = qJD(1) * t90;
t58 = t167 * t87 + t141;
t86 = t100 * t167 + t103 * t101;
t80 = t86 * qJD(1);
t111 = pkin(4) * t80 + t58;
t188 = qJD(4) + t111;
t187 = Ifges(5,4) - Ifges(4,5);
t186 = Ifges(5,5) - Ifges(4,6);
t77 = qJD(5) + t80;
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t136 = t167 * t101;
t130 = qJD(1) * t136;
t139 = t100 * t103;
t135 = qJD(3) * t139;
t72 = qJD(1) * t135 - qJD(3) * t130;
t115 = qJ(4) * t72 - qJD(4) * t80;
t171 = pkin(3) + pkin(7);
t82 = t86 * qJD(3);
t73 = qJD(1) * t82;
t14 = t171 * t73 + t115;
t109 = t86 * qJD(2);
t59 = -t103 * t87 + t167 * t88;
t30 = qJD(1) * t109 + qJD(3) * t59;
t16 = -t72 * pkin(4) + t30;
t137 = -pkin(2) * t101 - pkin(1);
t89 = qJD(1) * t137 + qJD(2);
t108 = -qJ(4) * t80 + t89;
t79 = qJD(1) * t139 - t130;
t23 = t171 * t79 + t108;
t27 = -qJD(3) * t171 + t188;
t5 = -t102 * t23 + t104 * t27;
t1 = qJD(5) * t5 + t102 * t16 + t104 * t14;
t6 = t102 * t27 + t104 * t23;
t2 = -qJD(5) * t6 - t102 * t14 + t104 * t16;
t185 = t1 * t102 + t104 * t2;
t184 = -t58 - qJD(4);
t63 = -qJD(3) * t102 + t104 * t79;
t40 = qJD(5) * t63 + t102 * t73;
t64 = qJD(3) * t104 + t102 * t79;
t41 = -qJD(5) * t64 + t104 * t73;
t183 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t40 + Ifges(6,6) * t41;
t182 = (m(3) * qJ(2) + mrSges(3,3)) * (t100 ^ 2 + t101 ^ 2);
t134 = qJD(3) * t167;
t94 = qJD(2) * t136;
t45 = (qJD(2) * t100 + qJD(3) * t91) * t103 + t90 * t134 - t94;
t138 = qJD(1) * qJD(2);
t133 = t100 * t138;
t152 = qJD(1) * t94 - t87 * t134;
t28 = qJD(3) * (-qJD(4) + t141) + t103 * t133 - t152;
t149 = Ifges(6,4) * t102;
t118 = Ifges(6,2) * t104 + t149;
t148 = Ifges(6,4) * t104;
t120 = Ifges(6,1) * t102 + t148;
t123 = mrSges(6,1) * t104 - mrSges(6,2) * t102;
t125 = t102 * t5 - t104 * t6;
t144 = Ifges(6,6) * t104;
t147 = Ifges(6,5) * t102;
t169 = -t102 / 0.2e1;
t174 = -t77 / 0.2e1;
t176 = -t64 / 0.2e1;
t178 = -t63 / 0.2e1;
t166 = Ifges(6,4) * t64;
t18 = Ifges(6,2) * t63 + Ifges(6,6) * t77 + t166;
t62 = Ifges(6,4) * t63;
t19 = Ifges(6,1) * t64 + Ifges(6,5) * t77 + t62;
t170 = t79 * pkin(4);
t54 = -qJD(3) * qJ(4) - t59;
t31 = -t54 - t170;
t181 = (t144 + t147) * t174 + t118 * t178 + t120 * t176 + t31 * t123 + t125 * mrSges(6,3) + t19 * t169 - t104 * t18 / 0.2e1;
t180 = t40 / 0.2e1;
t179 = t41 / 0.2e1;
t177 = t63 / 0.2e1;
t175 = t64 / 0.2e1;
t173 = t77 / 0.2e1;
t168 = t104 / 0.2e1;
t60 = t103 * t91 + t167 * t90;
t163 = t30 * t60;
t162 = t63 * Ifges(6,6);
t161 = t64 * Ifges(6,5);
t160 = t77 * Ifges(6,3);
t157 = -Ifges(4,4) - Ifges(5,6);
t65 = -qJD(3) * mrSges(4,2) - t79 * mrSges(4,3);
t67 = t79 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t155 = t65 - t67;
t154 = -qJD(3) * t158 - t189 * t80;
t33 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t153 = -t67 + t33;
t146 = Ifges(6,5) * t104;
t145 = Ifges(6,6) * t102;
t143 = qJ(4) * t79;
t142 = t102 * t82;
t140 = t104 * t82;
t128 = t1 * t104 - t2 * t102;
t126 = t6 * t102 + t5 * t104;
t61 = -t103 * t90 + t167 * t91;
t124 = -t60 * t72 - t61 * t73;
t122 = mrSges(6,1) * t102 + mrSges(6,2) * t104;
t121 = Ifges(6,1) * t104 - t149;
t119 = -Ifges(6,2) * t102 + t148;
t112 = -qJ(4) * t86 + t137;
t85 = -t136 + t139;
t37 = t171 * t85 + t112;
t47 = pkin(4) * t86 + t60;
t13 = t102 * t47 + t104 * t37;
t12 = -t102 * t37 + t104 * t47;
t43 = -mrSges(6,2) * t77 + mrSges(6,3) * t63;
t44 = mrSges(6,1) * t77 - mrSges(6,3) * t64;
t116 = -t102 * t43 - t104 * t44;
t81 = -t101 * t134 + t135;
t114 = qJ(4) * t81 - qJD(4) * t86;
t46 = qJD(3) * t61 + t109;
t75 = Ifges(4,4) * t79;
t74 = Ifges(5,6) * t79;
t71 = Ifges(6,3) * t72;
t70 = t72 * mrSges(4,2);
t69 = t72 * mrSges(5,3);
t57 = pkin(3) * t85 + t112;
t56 = -mrSges(5,2) * t79 - mrSges(5,3) * t80;
t55 = pkin(3) * t80 + t143;
t53 = -qJD(3) * pkin(3) - t184;
t52 = t80 * Ifges(4,1) + Ifges(4,5) * qJD(3) - t75;
t51 = t80 * Ifges(4,4) - t79 * Ifges(4,2) + Ifges(4,6) * qJD(3);
t50 = Ifges(5,4) * qJD(3) - t80 * Ifges(5,2) + t74;
t49 = Ifges(5,5) * qJD(3) - t80 * Ifges(5,6) + t79 * Ifges(5,3);
t48 = -t85 * pkin(4) + t61;
t42 = pkin(3) * t79 + t108;
t39 = t59 - t170;
t34 = pkin(3) * t82 + t114;
t32 = t171 * t80 + t143;
t29 = (-qJD(3) * t88 - t133) * t103 + t152;
t26 = pkin(3) * t73 + t115;
t25 = -t81 * pkin(4) + t46;
t24 = -pkin(4) * t82 - t45;
t22 = t171 * t82 + t114;
t21 = mrSges(6,2) * t72 + mrSges(6,3) * t41;
t20 = -mrSges(6,1) * t72 - mrSges(6,3) * t40;
t17 = t160 + t161 + t162;
t15 = -pkin(4) * t73 - t28;
t11 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t10 = t102 * t39 + t104 * t32;
t9 = -t102 * t32 + t104 * t39;
t8 = t40 * Ifges(6,1) + t41 * Ifges(6,4) - t72 * Ifges(6,5);
t7 = t40 * Ifges(6,4) + t41 * Ifges(6,2) - t72 * Ifges(6,6);
t4 = -qJD(5) * t13 - t102 * t22 + t104 * t25;
t3 = qJD(5) * t12 + t102 * t25 + t104 * t22;
t35 = [(-t71 / 0.2e1 - t26 * mrSges(5,3) + t157 * t73 + t189 * t30 + (-Ifges(6,3) / 0.2e1 - Ifges(4,1) - Ifges(5,2)) * t72 + t183) * t86 + 0.2e1 * t182 * t138 + m(6) * (t1 * t13 + t12 * t2 + t15 * t48 + t24 * t31 + t3 * t6 + t4 * t5) + t89 * (mrSges(4,1) * t82 - mrSges(4,2) * t81) - t80 * (Ifges(5,2) * t81 + Ifges(5,6) * t82) / 0.2e1 + t79 * (Ifges(5,6) * t81 + Ifges(5,3) * t82) / 0.2e1 + t42 * (-mrSges(5,2) * t82 + mrSges(5,3) * t81) + t82 * t49 / 0.2e1 + t80 * (-Ifges(4,1) * t81 - Ifges(4,4) * t82) / 0.2e1 - t79 * (-Ifges(4,4) * t81 - Ifges(4,2) * t82) / 0.2e1 - t82 * t51 / 0.2e1 + t81 * t50 / 0.2e1 + t57 * (-t73 * mrSges(5,2) + t69) + t34 * t56 + t3 * t43 + t4 * t44 + t48 * t11 + t24 * t33 + t12 * t20 + t13 * t21 - (t52 + t17) * t81 / 0.2e1 + (t186 * t82 + t187 * t81) * qJD(3) / 0.2e1 + (-t58 * t81 - t59 * t82 + t124) * mrSges(4,3) + (-t53 * t81 + t54 * t82 + t124) * mrSges(5,1) + (Ifges(6,5) * t142 + Ifges(6,6) * t140 - Ifges(6,3) * t81) * t173 + (Ifges(6,1) * t142 + Ifges(6,4) * t140 - Ifges(6,5) * t81) * t175 + (Ifges(6,4) * t142 + Ifges(6,2) * t140 - Ifges(6,6) * t81) * t177 + t137 * (t73 * mrSges(4,1) - t70) + t18 * t140 / 0.2e1 + t6 * (mrSges(6,2) * t81 + mrSges(6,3) * t140) + t19 * t142 / 0.2e1 + t5 * (-mrSges(6,1) * t81 - mrSges(6,3) * t142) + t31 * (-mrSges(6,1) * t140 + mrSges(6,2) * t142) - t154 * t46 - t155 * t45 + m(4) * (t29 * t61 - t45 * t59 + t46 * t58 + t163) + m(5) * (t26 * t57 - t28 * t61 + t34 * t42 + t45 * t54 + t46 * t53 + t163) + (t120 * t180 + t118 * t179 - t15 * t123 + t7 * t168 + t102 * t8 / 0.2e1 + t28 * mrSges(5,1) - t29 * mrSges(4,3) - t26 * mrSges(5,2) + (Ifges(4,2) + Ifges(5,3)) * t73 + (-t147 / 0.2e1 - t144 / 0.2e1 - t157) * t72 + t128 * mrSges(6,3) + ((-t145 + t146) * t173 + t119 * t177 + t121 * t175 + t31 * t122 + t19 * t168 + t18 * t169 - t126 * mrSges(6,3)) * qJD(5)) * t85; -t102 * t20 + t104 * t21 + t69 - t70 - t158 * t73 + t116 * qJD(5) + (t65 + t153) * t79 + (t116 + t154) * t80 - m(4) * (t58 * t80 - t59 * t79) - t182 * qJD(1) ^ 2 + (-t77 * t126 + t31 * t79 + t128) * m(6) + (-t53 * t80 - t54 * t79 + t26) * m(5); (qJ(4) * t15 - t10 * t6 + t188 * t31 - t5 * t9) * m(6) + t111 * t33 - t55 * t56 - t10 * t43 - t9 * t44 - t28 * mrSges(5,3) - t29 * mrSges(4,2) + qJ(4) * t11 - t185 * mrSges(6,3) + (-t54 * mrSges(5,1) + t59 * mrSges(4,3) - t89 * mrSges(4,1) + t51 / 0.2e1 + t42 * mrSges(5,2) - t49 / 0.2e1 + (Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t80 + (Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(3) + (-Ifges(4,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(5,3) / 0.2e1) * t79 + t181) * t80 + t181 * qJD(5) + t8 * t168 + (-pkin(3) * t30 - qJ(4) * t28 + t184 * t54 - t42 * t55 - t53 * t59) * m(5) - (t102 * t21 + t104 * t20 + m(6) * t185 + (-m(6) * t125 - t102 * t44 + t104 * t43) * qJD(5)) * t171 + (-mrSges(5,1) * qJ(4) + t186) * t73 + (pkin(3) * mrSges(5,1) - t146 / 0.2e1 + t145 / 0.2e1 + t187) * t72 + t15 * t122 + t7 * t169 + t119 * t179 + t121 * t180 + t153 * qJD(4) + t154 * t59 + t155 * t58 + t158 * t30 + (t53 * mrSges(5,1) + t58 * mrSges(4,3) + t52 / 0.2e1 + t17 / 0.2e1 - t50 / 0.2e1 - t75 / 0.2e1 - t74 / 0.2e1 - t42 * mrSges(5,3) - t6 * mrSges(6,2) + t5 * mrSges(6,1) + t161 / 0.2e1 + t160 / 0.2e1 + t162 / 0.2e1 + t89 * mrSges(4,2) + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * qJD(3)) * t79; -t72 * mrSges(5,1) + t80 * t56 - t153 * qJD(3) + (t43 * t77 + t20) * t104 + (-t44 * t77 + t21) * t102 + (-qJD(3) * t31 - t77 * t125 + t185) * m(6) + (qJD(3) * t54 + t42 * t80 + t30) * m(5); -t71 - t31 * (mrSges(6,1) * t64 + mrSges(6,2) * t63) + (Ifges(6,1) * t63 - t166) * t176 + t18 * t175 + (Ifges(6,5) * t63 - Ifges(6,6) * t64) * t174 - t5 * t43 + t6 * t44 + (t5 * t63 + t6 * t64) * mrSges(6,3) + (-Ifges(6,2) * t64 + t19 + t62) * t178 + t183;];
tauc = t35(:);
