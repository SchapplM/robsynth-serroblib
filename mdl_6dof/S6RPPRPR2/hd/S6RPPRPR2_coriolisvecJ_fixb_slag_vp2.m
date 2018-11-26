% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:39:59
% EndTime: 2018-11-23 15:40:02
% DurationCPUTime: 2.91s
% Computational Cost: add. (3568->373), mult. (8701->495), div. (0->0), fcn. (5952->8), ass. (0->171)
t175 = mrSges(6,2) - mrSges(5,1);
t174 = mrSges(5,3) + mrSges(6,1);
t107 = sin(pkin(10));
t109 = cos(pkin(10));
t146 = qJD(1) * (t107 ^ 2 + t109 ^ 2);
t206 = mrSges(4,3) * t146;
t112 = sin(qJ(4));
t160 = pkin(7) * qJD(1);
t101 = sin(pkin(9)) * pkin(1) + qJ(3);
t95 = t101 * qJD(1);
t78 = t107 * qJD(2) + t109 * t95;
t66 = t109 * t160 + t78;
t157 = t112 * t66;
t185 = cos(qJ(4));
t104 = t109 * qJD(2);
t65 = t104 + (-t95 - t160) * t107;
t36 = -t185 * t65 + t157;
t94 = t107 * t185 + t112 * t109;
t87 = t94 * qJD(1);
t123 = pkin(5) * t87 + t36;
t205 = t123 + qJD(5);
t204 = Ifges(6,4) - Ifges(5,5);
t203 = Ifges(6,5) - Ifges(5,6);
t202 = -qJD(5) - t36;
t84 = qJD(6) + t87;
t111 = sin(qJ(6));
t113 = cos(qJ(6));
t119 = t94 * qJD(3);
t37 = t112 * t65 + t185 * t66;
t20 = qJD(1) * t119 + qJD(4) * t37;
t151 = t185 * t109;
t145 = qJD(1) * t151;
t155 = t107 * t112;
t150 = qJD(4) * t155;
t79 = qJD(1) * t150 - qJD(4) * t145;
t15 = -t79 * pkin(5) + t20;
t127 = qJ(5) * t79 - qJD(5) * t87;
t189 = pkin(4) + pkin(8);
t89 = t94 * qJD(4);
t80 = qJD(1) * t89;
t18 = t189 * t80 + t127;
t16 = -qJD(4) * t189 + t205;
t124 = -cos(pkin(9)) * pkin(1) - pkin(3) * t109 - pkin(2);
t85 = qJD(1) * t124 + qJD(3);
t118 = -qJ(5) * t87 + t85;
t86 = qJD(1) * t155 - t145;
t31 = t189 * t86 + t118;
t5 = -t111 * t31 + t113 * t16;
t1 = qJD(6) * t5 + t111 * t15 + t113 * t18;
t6 = t111 * t16 + t113 * t31;
t2 = -qJD(6) * t6 - t111 * t18 + t113 * t15;
t201 = t1 * t111 + t113 * t2;
t67 = -qJD(4) * t111 + t113 * t86;
t49 = qJD(6) * t67 + t111 * t80;
t68 = qJD(4) * t113 + t111 * t86;
t50 = -qJD(6) * t68 + t113 * t80;
t200 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t49 + Ifges(7,6) * t50;
t149 = qJD(4) * t185;
t153 = qJD(3) * t107;
t172 = pkin(7) + t101;
t90 = t172 * t107;
t91 = t172 * t109;
t98 = qJD(3) * t151;
t38 = t112 * (qJD(4) * t91 + t153) + t90 * t149 - t98;
t148 = qJD(1) * t153;
t171 = qJD(1) * t98 + t65 * t149;
t17 = qJD(4) * (-qJD(5) + t157) + t112 * t148 - t171;
t166 = Ifges(7,4) * t111;
t133 = Ifges(7,2) * t113 + t166;
t165 = Ifges(7,4) * t113;
t135 = Ifges(7,1) * t111 + t165;
t138 = mrSges(7,1) * t113 - mrSges(7,2) * t111;
t140 = t111 * t5 - t113 * t6;
t161 = Ifges(7,6) * t113;
t164 = Ifges(7,5) * t111;
t187 = -t111 / 0.2e1;
t192 = -t84 / 0.2e1;
t194 = -t68 / 0.2e1;
t196 = -t67 / 0.2e1;
t188 = t86 * pkin(5);
t34 = -qJD(4) * qJ(5) - t37;
t21 = -t34 - t188;
t184 = Ifges(7,4) * t68;
t27 = t67 * Ifges(7,2) + t84 * Ifges(7,6) + t184;
t64 = Ifges(7,4) * t67;
t28 = t68 * Ifges(7,1) + t84 * Ifges(7,5) + t64;
t199 = (t161 + t164) * t192 + t133 * t196 + t135 * t194 + t21 * t138 + t140 * mrSges(7,3) + t28 * t187 - t113 * t27 / 0.2e1;
t198 = t49 / 0.2e1;
t197 = t50 / 0.2e1;
t195 = t67 / 0.2e1;
t193 = t68 / 0.2e1;
t191 = t84 / 0.2e1;
t186 = t113 / 0.2e1;
t58 = t112 * t91 + t185 * t90;
t181 = t20 * t58;
t93 = -t151 + t155;
t180 = t20 * t93;
t179 = t67 * Ifges(7,6);
t178 = t68 * Ifges(7,5);
t177 = t84 * Ifges(7,3);
t173 = -Ifges(5,4) - Ifges(6,6);
t69 = -qJD(4) * mrSges(5,2) - t86 * mrSges(5,3);
t71 = mrSges(6,1) * t86 - qJD(4) * mrSges(6,3);
t170 = t69 - t71;
t169 = -qJD(4) * t175 - t174 * t87;
t42 = -mrSges(7,1) * t67 + mrSges(7,2) * t68;
t168 = -t71 + t42;
t163 = Ifges(7,5) * t113;
t162 = Ifges(7,6) * t111;
t159 = qJ(5) * t86;
t158 = t111 * t89;
t156 = t113 * t89;
t152 = t69 + t168;
t143 = t1 * t113 - t2 * t111;
t141 = t6 * t111 + t5 * t113;
t59 = -t112 * t90 + t185 * t91;
t139 = -t58 * t79 - t59 * t80;
t137 = mrSges(7,1) * t111 + mrSges(7,2) * t113;
t136 = Ifges(7,1) * t113 - t166;
t134 = -Ifges(7,2) * t111 + t165;
t131 = -t107 * (-t107 * t95 + t104) + t109 * t78;
t29 = -mrSges(7,1) * t79 - mrSges(7,3) * t49;
t30 = mrSges(7,2) * t79 + mrSges(7,3) * t50;
t130 = t111 * t30 + t113 * t29;
t120 = -qJ(5) * t94 + t124;
t40 = t189 * t93 + t120;
t43 = pkin(5) * t94 + t58;
t12 = t111 * t43 + t113 * t40;
t11 = -t111 * t40 + t113 * t43;
t51 = -mrSges(7,2) * t84 + mrSges(7,3) * t67;
t52 = mrSges(7,1) * t84 - mrSges(7,3) * t68;
t129 = -t111 * t52 + t113 * t51;
t128 = -t111 * t51 - t113 * t52;
t88 = -t109 * t149 + t150;
t126 = qJ(5) * t88 - qJD(5) * t94;
t121 = -t128 - t169;
t117 = -qJD(6) * t140 + t201;
t39 = qJD(4) * t59 + t119;
t83 = Ifges(5,4) * t86;
t82 = Ifges(6,6) * t86;
t75 = Ifges(7,3) * t79;
t74 = t79 * mrSges(5,2);
t73 = t79 * mrSges(6,3);
t61 = -mrSges(6,2) * t86 - mrSges(6,3) * t87;
t60 = pkin(4) * t87 + t159;
t57 = t87 * Ifges(5,1) + Ifges(5,5) * qJD(4) - t83;
t56 = t87 * Ifges(5,4) - t86 * Ifges(5,2) + Ifges(5,6) * qJD(4);
t55 = Ifges(6,4) * qJD(4) - t87 * Ifges(6,2) + t82;
t54 = Ifges(6,5) * qJD(4) - t87 * Ifges(6,6) + t86 * Ifges(6,3);
t53 = pkin(4) * t93 + t120;
t46 = pkin(4) * t89 + t126;
t45 = pkin(4) * t86 + t118;
t44 = -t93 * pkin(5) + t59;
t41 = t189 * t87 + t159;
t35 = pkin(4) * t80 + t127;
t33 = -qJD(4) * pkin(4) - t202;
t32 = t189 * t89 + t126;
t26 = t177 + t178 + t179;
t25 = -t88 * pkin(5) + t39;
t24 = -pkin(5) * t89 - t38;
t23 = t37 - t188;
t19 = (-qJD(4) * t66 - t148) * t112 + t171;
t14 = -pkin(5) * t80 - t17;
t13 = -mrSges(7,1) * t50 + mrSges(7,2) * t49;
t10 = t49 * Ifges(7,1) + t50 * Ifges(7,4) - t79 * Ifges(7,5);
t9 = t49 * Ifges(7,4) + t50 * Ifges(7,2) - t79 * Ifges(7,6);
t8 = t111 * t23 + t113 * t41;
t7 = -t111 * t41 + t113 * t23;
t4 = -qJD(6) * t12 - t111 * t32 + t113 * t25;
t3 = qJD(6) * t11 + t111 * t25 + t113 * t32;
t22 = [(m(4) * (t101 * t146 + t131) + 0.2e1 * t206) * qJD(3) - (t57 + t26) * t88 / 0.2e1 + (-t36 * t88 - t37 * t89 + t139) * mrSges(5,3) + (-t33 * t88 + t34 * t89 + t139) * mrSges(6,1) + t124 * (t80 * mrSges(5,1) - t74) + t53 * (-t80 * mrSges(6,2) + t73) + t46 * t61 + t3 * t51 + t4 * t52 + t24 * t42 + t44 * t13 + t11 * t29 + t12 * t30 + m(7) * (t1 * t12 + t11 * t2 + t14 * t44 + t21 * t24 + t3 * t6 + t4 * t5) + (Ifges(7,5) * t158 + Ifges(7,6) * t156 - Ifges(7,3) * t88) * t191 + (Ifges(7,1) * t158 + Ifges(7,4) * t156 - Ifges(7,5) * t88) * t193 + t27 * t156 / 0.2e1 + t6 * (mrSges(7,2) * t88 + mrSges(7,3) * t156) + t28 * t158 / 0.2e1 + t5 * (-mrSges(7,1) * t88 - mrSges(7,3) * t158) + t21 * (-mrSges(7,1) * t156 + mrSges(7,2) * t158) - t169 * t39 - t170 * t38 + m(5) * (t19 * t59 + t36 * t39 - t37 * t38 + t181) + m(6) * (-t17 * t59 + t33 * t39 + t34 * t38 + t35 * t53 + t45 * t46 + t181) + (t17 * mrSges(6,1) - t19 * mrSges(5,3) - t35 * mrSges(6,2) - t14 * t138 + t135 * t198 + t133 * t197 + t111 * t10 / 0.2e1 + t9 * t186 + (Ifges(5,2) + Ifges(6,3)) * t80 + (-t164 / 0.2e1 - t161 / 0.2e1 - t173) * t79 + t143 * mrSges(7,3) + (t27 * t187 + t28 * t186 + t21 * t137 + (-t162 + t163) * t191 + t136 * t193 + t134 * t195 - t141 * mrSges(7,3)) * qJD(6)) * t93 + (Ifges(7,4) * t158 + Ifges(7,2) * t156 - Ifges(7,6) * t88) * t195 + t88 * t55 / 0.2e1 + t89 * t54 / 0.2e1 + t86 * (Ifges(6,6) * t88 + Ifges(6,3) * t89) / 0.2e1 + t45 * (-mrSges(6,2) * t89 + mrSges(6,3) * t88) - t87 * (Ifges(6,2) * t88 + Ifges(6,6) * t89) / 0.2e1 - t89 * t56 / 0.2e1 + t85 * (mrSges(5,1) * t89 - mrSges(5,2) * t88) - t86 * (-Ifges(5,4) * t88 - Ifges(5,2) * t89) / 0.2e1 + t87 * (-Ifges(5,1) * t88 - Ifges(5,4) * t89) / 0.2e1 + (-t75 / 0.2e1 - t35 * mrSges(6,3) + t173 * t80 + t174 * t20 + (-Ifges(5,1) - Ifges(6,2) - Ifges(7,3) / 0.2e1) * t79 + t200) * t94 + (t203 * t89 + t204 * t88) * qJD(4) / 0.2e1; (-t174 * t80 + t13) * t94 - t152 * t88 + t121 * t89 + (t129 * qJD(6) - t174 * t79 + t130) * t93 + m(5) * (t19 * t94 + t36 * t89 - t37 * t88 + t180) + m(6) * (-t17 * t94 + t33 * t89 + t34 * t88 + t180) + m(7) * (t117 * t93 + t14 * t94 + t141 * t89 - t21 * t88); -t111 * t29 + t113 * t30 + t73 - t74 - t175 * t80 + t128 * qJD(6) + t152 * t86 - t121 * t87 - m(5) * (t36 * t87 - t37 * t86) + (-m(4) * t131 - t206) * qJD(1) + (-t141 * t84 + t21 * t86 + t143) * m(7) + (-t33 * t87 - t34 * t86 + t35) * m(6); (qJ(5) * t14 + t205 * t21 - t5 * t7 - t6 * t8) * m(7) - t201 * mrSges(7,3) + t123 * t42 + t14 * t137 - t60 * t61 - t8 * t51 - t7 * t52 - t17 * mrSges(6,3) - t19 * mrSges(5,2) + qJ(5) * t13 + t10 * t186 + t9 * t187 + t168 * qJD(5) + t169 * t37 + t170 * t36 + t175 * t20 + (-t55 / 0.2e1 + t57 / 0.2e1 + t26 / 0.2e1 + t177 / 0.2e1 + t178 / 0.2e1 + t179 / 0.2e1 - t82 / 0.2e1 - t83 / 0.2e1 - t45 * mrSges(6,3) + t36 * mrSges(5,3) + t33 * mrSges(6,1) + t85 * mrSges(5,2) + t5 * mrSges(7,1) - t6 * mrSges(7,2) + (Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * qJD(4)) * t86 + (t37 * mrSges(5,3) - t34 * mrSges(6,1) + t45 * mrSges(6,2) - t54 / 0.2e1 - t85 * mrSges(5,1) + t56 / 0.2e1 + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t87 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t86 + t199) * t87 + t199 * qJD(6) + t134 * t197 + t136 * t198 - (t130 + m(7) * t201 + (-m(7) * t140 + t129) * qJD(6)) * t189 + (-pkin(4) * t20 - qJ(5) * t17 + t202 * t34 - t33 * t37 - t45 * t60) * m(6) + (-mrSges(6,1) * qJ(5) + t203) * t80 + (pkin(4) * mrSges(6,1) - t163 / 0.2e1 + t162 / 0.2e1 + t204) * t79; -t79 * mrSges(6,1) + t87 * t61 - t168 * qJD(4) + (t51 * t84 + t29) * t113 + (-t52 * t84 + t30) * t111 + (-qJD(4) * t21 - t140 * t87 + t117) * m(7) + (qJD(4) * t34 + t45 * t87 + t20) * m(6); -t75 - t21 * (mrSges(7,1) * t68 + mrSges(7,2) * t67) + (Ifges(7,1) * t67 - t184) * t194 + t27 * t193 + (Ifges(7,5) * t67 - Ifges(7,6) * t68) * t192 - t5 * t51 + t6 * t52 + (t5 * t67 + t6 * t68) * mrSges(7,3) + (-Ifges(7,2) * t68 + t28 + t64) * t196 + t200;];
tauc  = t22(:);
