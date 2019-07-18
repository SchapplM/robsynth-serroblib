% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:49
% EndTime: 2019-07-18 13:28:55
% DurationCPUTime: 1.83s
% Computational Cost: add. (1791->260), mult. (5380->393), div. (0->0), fcn. (4314->8), ass. (0->154)
t79 = sin(qJ(2));
t144 = qJD(1) * t79;
t77 = sin(qJ(4));
t81 = cos(qJ(3));
t154 = t77 * t81;
t182 = cos(qJ(4));
t78 = sin(qJ(3));
t126 = t78 * t144;
t93 = -qJD(3) * pkin(2) + t126;
t41 = t144 * t154 + t182 * t93;
t189 = t41 * qJD(4);
t133 = qJD(3) + qJD(4);
t59 = t182 * t78 + t154;
t188 = t133 * t59;
t27 = t188 * qJD(2);
t122 = t182 * t81;
t108 = qJD(1) * t122;
t99 = t79 * t108;
t42 = -t77 * t93 + t99;
t82 = cos(qJ(2));
t143 = qJD(1) * t82;
t63 = -pkin(2) * qJD(2) * t81 - t143;
t76 = sin(qJ(5));
t80 = cos(qJ(5));
t100 = t42 * t76 - t63 * t80;
t187 = t100 * qJD(5);
t134 = qJD(2) * qJD(3);
t186 = -0.2e1 * t134;
t55 = t59 * qJD(2);
t40 = t76 * t133 + t80 * t55;
t112 = qJD(3) * t126;
t140 = qJD(2) * t82;
t151 = t79 * t81;
t48 = (-qJD(3) * t151 - t78 * t140) * qJD(1);
t87 = -t182 * t112 + t77 * t48 - t189;
t98 = t82 * t108;
t12 = qJD(2) * t98 + t87;
t24 = t42 * t80 + t63 * t76;
t138 = qJD(5) * t24;
t139 = qJD(3) * t78;
t57 = (pkin(2) * t139 + t144) * qJD(2);
t8 = -t12 * t76 + t80 * t57 - t138;
t185 = t40 * t41 - t8;
t74 = t78 ^ 2;
t75 = t81 ^ 2;
t146 = t74 + t75;
t184 = (-0.1e1 + t146) * t79;
t117 = t182 * qJD(4);
t183 = t182 * qJD(3) + t117;
t155 = t77 * t78;
t106 = t133 * t155;
t119 = qJD(2) * t182;
t69 = t81 * t119;
t148 = t133 * t69;
t26 = qJD(2) * t106 - t148;
t15 = t40 * qJD(5) - t76 * t26;
t7 = t12 * t80 + t57 * t76 - t187;
t4 = t7 * t80;
t115 = t80 * t133;
t137 = qJD(5) * t76;
t14 = -qJD(5) * t115 + t55 * t137 + t80 * t26;
t181 = t14 * t76;
t180 = t15 * t80;
t94 = t122 - t155;
t179 = t27 * t94;
t178 = t27 * t76;
t177 = t27 * t80;
t176 = t27 * t81;
t38 = t55 * t76 - t115;
t175 = t38 * t41;
t142 = qJD(2) * t78;
t127 = t77 * t142;
t53 = -t69 + t127;
t52 = qJD(5) + t53;
t174 = t38 * t52;
t173 = t38 * t76;
t172 = t38 * t80;
t171 = t40 * t38;
t169 = t40 * t52;
t168 = t40 * t76;
t167 = t40 * t80;
t44 = t94 * t144;
t166 = t41 * t44;
t45 = t59 * t143;
t165 = t41 * t45;
t164 = t52 * t55;
t163 = t53 * t76;
t162 = t53 * t80;
t161 = t55 * t53;
t160 = t59 * t76;
t159 = t59 * t80;
t158 = t63 * t53;
t157 = t63 * t55;
t153 = t78 * t79;
t152 = t78 * t81;
t84 = qJD(2) ^ 2;
t150 = t84 * t82;
t149 = -t77 * t112 - t182 * t48;
t147 = t74 - t75;
t83 = qJD(3) ^ 2;
t145 = t83 + t84;
t141 = qJD(2) * t79;
t136 = qJD(5) * t80;
t135 = qJD(5) * t81;
t131 = t100 * t162 - t24 * t163 + t4;
t129 = t84 * t152;
t128 = pkin(2) * t142;
t125 = t76 * t182;
t124 = t80 * t182;
t120 = qJD(1) * t140;
t110 = t81 * t120;
t92 = t77 * t110 + t149;
t13 = t42 * qJD(4) + t92;
t123 = t182 * t13;
t121 = t41 * (-t52 + t53);
t116 = t52 * t80;
t113 = t82 * t186;
t111 = t134 * t152;
t109 = t13 * t76 + t41 * t136 + t24 * t55;
t107 = t52 * t117;
t19 = t78 * t82 * t119 + (t77 * t140 + t183 * t79) * t81 - t133 * t77 * t153;
t49 = t59 * t79;
t104 = t13 * t49 + t19 * t41;
t103 = -t100 * t80 + t24 * t76;
t102 = -t100 * t76 - t24 * t80;
t101 = t167 + t173;
t50 = t94 * t79;
t37 = t50 * t80 - t76 * t82;
t36 = -t50 * t76 - t80 * t82;
t97 = t100 * t55 - t13 * t80 + t41 * t137;
t32 = -t183 * t81 + t106;
t96 = t59 * t136 - t32 * t76;
t95 = -t59 * t137 - t32 * t80;
t91 = t82 * t94;
t89 = -t103 * qJD(5) - t8 * t76;
t88 = t102 * qJD(5) - t7 * t76 - t8 * t80;
t47 = qJD(1) * t91;
t46 = t59 * t144;
t31 = t76 * t144 + t47 * t80;
t30 = t80 * t144 - t47 * t76;
t29 = t76 * t128 - t46 * t80;
t28 = t80 * t128 + t46 * t76;
t22 = -t53 ^ 2 + t55 ^ 2;
t21 = t55 * t133 - t27;
t20 = t148 + (-t127 + t53) * t133;
t18 = qJD(2) * t91 - t188 * t79;
t10 = -t37 * qJD(5) + t80 * t141 - t18 * t76;
t9 = t36 * qJD(5) + t76 * t141 + t18 * t80;
t6 = t52 * t116 - t40 * t55 + t178;
t5 = -t52 ^ 2 * t76 + t38 * t55 + t177;
t3 = t52 * t173 - t180;
t2 = t40 * t116 - t181;
t1 = (-t14 - t174) * t80 + (-t15 - t169) * t76;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t79, -t150, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t113 - t145 * t151, t81 * t113 + t145 * t153, t146 * t150, 0.2e1 * t120 * t184, 0, 0, 0, 0, 0, 0, -t19 * t133 + t53 * t141 - t82 * t27, -t18 * t133 + t55 * t141 + t82 * t26, -t18 * t53 + t19 * t55 - t26 * t49 - t27 * t50, t12 * t50 + t63 * t141 + t18 * t42 - t57 * t82 + t104, 0, 0, 0, 0, 0, 0, t10 * t52 + t15 * t49 + t19 * t38 + t27 * t36, -t14 * t49 + t19 * t40 - t27 * t37 - t52 * t9, -t10 * t40 + t14 * t36 - t15 * t37 - t38 * t9, -t10 * t100 + t24 * t9 + t36 * t8 + t37 * t7 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t111, t147 * t186, t83 * t81, -0.2e1 * t111, -t83 * t78, 0, 0, 0, 0, -qJD(1) ^ 2 * t82 * t184, -t26 * t59 - t32 * t55, -t188 * t55 - t26 * t94 - t27 * t59 + t32 * t53, -t32 * t133, t188 * t53 - t179, -t188 * t133, 0, -t57 * t94 + t63 * t188 - t53 * t144 + t45 * t133 + (t53 * t139 - t176) * pkin(2), t57 * t59 - t63 * t32 - t55 * t144 + t47 * t133 + (t55 * t139 + t26 * t81) * pkin(2), t12 * t94 + t13 * t59 - t188 * t42 - t32 * t41 - t45 * t55 + t47 * t53, -t63 * t144 - t165 - t42 * t47 + (t63 * t139 - t57 * t81) * pkin(2), -t14 * t159 + t40 * t95, (t168 + t172) * t32 + (t181 - t180 + (-t167 + t173) * qJD(5)) * t59, t14 * t94 + t27 * t159 + t188 * t40 + t52 * t95, t15 * t160 + t38 * t96, t15 * t94 - t27 * t160 - t188 * t38 - t52 * t96, t188 * t52 - t179, t13 * t160 - t100 * t188 - t30 * t52 - t38 * t45 - t94 * t8 + t96 * t41 + (-t80 * t176 + (t76 * t135 + t80 * t139) * t52) * pkin(2), t13 * t159 - t24 * t188 + t31 * t52 - t40 * t45 + t94 * t7 + t95 * t41 + (t76 * t176 + (t80 * t135 - t76 * t139) * t52) * pkin(2), t30 * t40 + t31 * t38 + t103 * t32 + t88 * t59 + (-t101 * t139 + (-t14 * t80 + t15 * t76 + (-t168 + t172) * qJD(5)) * t81) * pkin(2), t100 * t30 - t24 * t31 - t165 + (t103 * t139 + t81 * t88) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t147 * t84, 0, t129, 0, 0, 0, 0, 0, 0, t161, t22, t20, -t161, t21, 0, -qJD(4) * t99 - t53 * t128 - t157 + t44 * t133 + (-t110 + (t126 + (-0.2e1 * qJD(3) - qJD(4)) * pkin(2)) * qJD(4)) * t77 - t149, t158 + (-t78 * pkin(2) * t55 - t98) * qJD(2) - t87 + (-pkin(2) * t117 - t46) * t133, (t42 - t44) * t55 + (t41 - t46) * t53 + (t182 * t26 - t27 * t77 + (-t182 * t53 + t55 * t77) * qJD(4)) * pkin(2), -t166 + t42 * t46 + (-t63 * t142 - t123 + t12 * t77 + (t182 * t42 + t41 * t77) * qJD(4)) * pkin(2), t2, t1, t6, t3, t5, -t164, t41 * t163 - t28 * t52 - t44 * t38 + (-t76 * t107 - t182 * t15 + (qJD(4) * t38 - t52 * t136 - t178) * t77) * pkin(2) + t97, t41 * t162 + t29 * t52 - t44 * t40 + (-t80 * t107 + t182 * t14 + (qJD(4) * t40 + t52 * t137 - t177) * t77) * pkin(2) + t109, t28 * t40 + t29 * t38 + ((-t124 * t38 + t125 * t40) * qJD(4) + (qJD(5) * t101 - t180 - t181) * t77) * pkin(2) + t89 + t131, t100 * t28 - t24 * t29 - t166 + (-t123 + (t100 * t125 + t124 * t24) * qJD(4) + (t4 + t89 + t189) * t77) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, t22, t20, -t161, t21, 0, t42 * qJD(3) - t157 - t92, -t41 * t133 - t12 + t158, 0, 0, t2, t1, t6, t3, t5, -t164, t121 * t76 - t38 * t42 + t97, t121 * t80 - t40 * t42 + t109, (-t175 + t187) * t80 + (-t138 + t185) * t76 + t131, (-t102 - t42) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, -t38 ^ 2 + t40 ^ 2, -t14 + t174, -t171, -t15 + t169, t27, t24 * t52 - t185, -t100 * t52 + t175 - t7, 0, 0;];
tauc_reg  = t11;
