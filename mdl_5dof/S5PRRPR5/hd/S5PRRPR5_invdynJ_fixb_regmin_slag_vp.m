% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:42
% EndTime: 2021-01-15 16:04:53
% DurationCPUTime: 2.42s
% Computational Cost: add. (2003->318), mult. (4901->457), div. (0->0), fcn. (3948->14), ass. (0->173)
t107 = sin(pkin(10));
t113 = sin(qJ(3));
t186 = qJD(2) * t113;
t116 = cos(qJ(3));
t200 = cos(pkin(10));
t157 = t200 * t116;
t90 = qJD(2) * t157;
t71 = t107 * t186 - t90;
t65 = qJD(5) + t71;
t221 = t65 - qJD(5);
t111 = qJ(4) + pkin(7);
t104 = qJ(3) + pkin(10);
t100 = cos(t104);
t110 = cos(pkin(5));
t117 = cos(qJ(2));
t201 = cos(pkin(9));
t159 = t201 * t117;
t108 = sin(pkin(9));
t114 = sin(qJ(2));
t195 = t108 * t114;
t68 = -t110 * t159 + t195;
t160 = t201 * t114;
t194 = t108 * t117;
t70 = t110 * t194 + t160;
t149 = g(1) * t70 + g(2) * t68;
t109 = sin(pkin(5));
t191 = t109 * t117;
t132 = g(3) * t191 - t149;
t126 = t132 * t100;
t179 = t113 * qJDD(2);
t146 = -qJDD(2) * t157 + t107 * t179;
t158 = t200 * t113;
t80 = t107 * t116 + t158;
t73 = t80 * qJD(3);
t38 = qJD(2) * t73 + t146;
t34 = qJDD(5) + t38;
t197 = t107 * t113;
t136 = t157 - t197;
t98 = pkin(3) * t116 + pkin(2);
t36 = -pkin(4) * t136 - pkin(8) * t80 - t98;
t220 = t36 * t34 - t126;
t74 = t80 * qJD(2);
t188 = qJD(1) * t110;
t182 = qJD(1) * qJD(2);
t59 = qJDD(2) * pkin(7) + (qJDD(1) * t114 + t117 * t182) * t109;
t123 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t188 + t59;
t187 = qJD(1) * t114;
t171 = t109 * t187;
t152 = t111 * qJD(2) + t171;
t144 = t152 * qJD(3);
t180 = t110 * qJDD(1);
t89 = t116 * t180;
t12 = qJDD(3) * pkin(3) - t123 * t113 - t116 * t144 + t89;
t13 = (-t144 + t180) * t113 + t123 * t116;
t3 = -t107 * t13 + t200 * t12;
t1 = -qJDD(3) * pkin(4) - t3;
t67 = t110 * t195 - t159;
t69 = t110 * t160 + t194;
t150 = g(1) * t67 - g(2) * t69;
t193 = t109 * t114;
t133 = -g(3) * t193 + t150;
t50 = t113 * t188 + t152 * t116;
t204 = t107 * t50;
t205 = qJD(3) * pkin(3);
t49 = -t152 * t113 + t116 * t188;
t47 = t49 + t205;
t18 = t200 * t47 - t204;
t14 = -qJD(3) * pkin(4) - t18;
t4 = t107 * t12 + t200 * t13;
t2 = qJDD(3) * pkin(8) + t4;
t162 = qJD(3) * t111;
t131 = -t113 * qJD(4) - t116 * t162;
t170 = qJD(1) * t191;
t66 = t116 * qJD(4) - t113 * t162;
t207 = t107 * t131 - t136 * t170 + t200 * t66;
t64 = -t98 * qJD(2) + qJD(4) - t170;
t24 = pkin(4) * t71 - pkin(8) * t74 + t64;
t84 = t111 * t116;
t54 = -t111 * t197 + t200 * t84;
t76 = t136 * qJD(3);
t219 = (qJD(5) * t24 + t2) * t136 + t1 * t80 + t14 * t76 + (-qJD(5) * t36 - t207) * t65 - t54 * t34 + t133;
t161 = t109 * t201;
t196 = t108 * t109;
t99 = sin(t104);
t137 = g(1) * (t100 * t196 + t67 * t99) + g(2) * (-t100 * t161 - t69 * t99) + g(3) * (t100 * t110 - t99 * t193);
t177 = pkin(3) * t186;
t215 = pkin(3) * t107;
t95 = pkin(8) + t215;
t218 = (pkin(4) * t74 + pkin(8) * t71 + qJD(5) * t95 + t177) * t65 + t1 + t137;
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t181 = qJD(2) * qJD(3);
t166 = t113 * t181;
t125 = t80 * qJDD(2) - t107 * t166;
t39 = qJD(3) * t90 + t125;
t57 = qJD(3) * t112 + t115 * t74;
t17 = t57 * qJD(5) - t115 * qJDD(3) + t112 * t39;
t118 = qJD(3) ^ 2;
t164 = qJDD(1) * t191;
t167 = t114 * t182;
t87 = t109 * t167;
t217 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t118 + (-g(3) * t117 + t167) * t109 + t149 + t164 - t87;
t77 = t110 * t116 - t113 * t193;
t216 = g(3) * t77;
t214 = pkin(3) * t113;
t183 = t115 * qJD(3);
t55 = t112 * t74 - t183;
t212 = t55 * t65;
t211 = t55 * t74;
t210 = t57 * t65;
t209 = t57 * t74;
t208 = t107 * t66 - t200 * t131 - t80 * t170;
t45 = t200 * t50;
t19 = t107 * t47 + t45;
t206 = qJD(2) * pkin(2);
t203 = t112 * t34;
t153 = t115 * t65;
t184 = qJD(5) * t112;
t16 = qJD(5) * t183 + t112 * qJDD(3) + t115 * t39 - t74 * t184;
t202 = t16 * t112;
t199 = qJD(5) * t80;
t192 = t109 * t116;
t190 = qJDD(1) - g(3);
t105 = t113 ^ 2;
t189 = -t116 ^ 2 + t105;
t185 = qJD(2) * t114;
t178 = t116 * qJDD(2);
t176 = t113 * t205;
t175 = t80 * t184;
t174 = t112 * t191;
t173 = t115 * t191;
t172 = t200 * pkin(3);
t169 = t109 * t185;
t168 = qJD(2) * t191;
t165 = t116 * t181;
t154 = t109 * t190;
t151 = t113 * t168;
t148 = pkin(4) * t73 - pkin(8) * t76 - t171 + t176;
t147 = t34 * t80 + t65 * t76;
t15 = qJD(3) * pkin(8) + t19;
t6 = t112 * t24 + t115 * t15;
t145 = t112 * t15 - t115 * t24;
t143 = t115 * t34 + (-t112 * t71 - t184) * t65;
t119 = qJD(2) ^ 2;
t142 = qJDD(2) * t117 - t114 * t119;
t43 = -t100 * t67 + t99 * t196;
t141 = -g(1) * t108 + t201 * g(2);
t78 = t110 * t113 + t114 * t192;
t30 = t107 * t77 + t200 * t78;
t140 = -t112 * t30 - t173;
t139 = -t115 * t30 + t174;
t135 = t142 * t109;
t83 = -t170 - t206;
t130 = -qJD(2) * t83 - t150 - t59;
t129 = t78 * qJD(3);
t128 = pkin(3) * t166 - t98 * qJDD(2) + qJDD(4) + t87;
t23 = t200 * t49 - t204;
t127 = -t95 * t34 + (t14 + t23) * t65;
t122 = -pkin(7) * qJDD(3) + (t170 + t83 - t206) * qJD(3);
t37 = t128 - t164;
t121 = -t129 - t151;
t96 = -t172 - pkin(4);
t62 = t100 * t193 + t110 * t99;
t53 = t107 * t84 + t111 * t158;
t48 = t77 * qJD(3) + t116 * t168;
t41 = t69 * t100 - t99 * t161;
t29 = t107 * t78 - t200 * t77;
t22 = t107 * t121 + t200 * t48;
t21 = t107 * t49 + t45;
t20 = t107 * t48 - t200 * t121;
t8 = pkin(4) * t38 - pkin(8) * t39 + t37;
t7 = t115 * t8;
t5 = [t190, 0, t135, (-qJDD(2) * t114 - t117 * t119) * t109, 0, 0, 0, 0, 0, t77 * qJDD(3) + t116 * t135 + (-t129 - 0.2e1 * t151) * qJD(3), -qJD(3) * t48 - qJDD(3) * t78 + (-t142 * t113 - t117 * t165) * t109, -qJD(3) * t20 - qJDD(3) * t29 + (-t117 * t38 + t71 * t185) * t109, -qJD(3) * t22 - qJDD(3) * t30 + (-t117 * t39 + t74 * t185) * t109, t20 * t74 - t22 * t71 + t29 * t39 - t30 * t38, -t18 * t20 + t19 * t22 - t29 * t3 + t30 * t4 - g(3) + (-t117 * t37 + t64 * t185) * t109, 0, 0, 0, 0, 0, (qJD(5) * t139 - t112 * t22 + t115 * t169) * t65 + t140 * t34 + t20 * t55 + t29 * t17, -(qJD(5) * t140 + t112 * t169 + t115 * t22) * t65 + t139 * t34 + t20 * t57 + t29 * t16; 0, qJDD(2), t190 * t191 + t149, -t114 * t154 - t150, qJDD(2) * t105 + 0.2e1 * t113 * t165, 0.2e1 * t113 * t178 - 0.2e1 * t189 * t181, qJDD(3) * t113 + t116 * t118, qJDD(3) * t116 - t113 * t118, 0, t122 * t113 + t116 * t217, -t113 * t217 + t122 * t116, -t71 * t171 - qJDD(3) * t53 - t37 * t136 - t38 * t98 + t64 * t73 - t126 + (t71 * t214 - t208) * qJD(3), -t74 * t171 - qJDD(3) * t54 + t37 * t80 - t39 * t98 + t64 * t76 + t132 * t99 + (t74 * t214 - t207) * qJD(3), t136 * t4 - t18 * t76 - t19 * t73 - t207 * t71 + t208 * t74 - t3 * t80 - t38 * t54 + t39 * t53 + t133, t4 * t54 - t3 * t53 - t37 * t98 + t64 * t176 - g(1) * (-t111 * t67 - t70 * t98) - g(2) * (t111 * t69 - t68 * t98) + t207 * t19 - t208 * t18 + (-t64 * t187 - g(3) * (t111 * t114 + t117 * t98)) * t109, -t57 * t175 + (t16 * t80 + t57 * t76) * t115, (-t112 * t57 - t115 * t55) * t76 + (-t202 - t115 * t17 + (t112 * t55 - t115 * t57) * qJD(5)) * t80, t115 * t147 - t136 * t16 - t175 * t65 + t57 * t73, -t112 * t147 + t136 * t17 - t199 * t153 - t55 * t73, -t136 * t34 + t65 * t73, t53 * t17 - t145 * t73 - t7 * t136 + t208 * t55 + (t148 * t65 + (t136 * t15 + t14 * t80 - t54 * t65) * qJD(5) + t220) * t115 + t219 * t112, t53 * t16 - t6 * t73 + t208 * t57 + ((-qJD(5) * t15 + t8) * t136 - t14 * t199 + (qJD(5) * t54 - t148) * t65 - t220) * t112 + t219 * t115; 0, 0, 0, 0, -t113 * t119 * t116, t189 * t119, t179, t178, qJDD(3), t130 * t113 + t141 * t192 - t216 + t89, g(3) * t78 + (-t141 * t109 - t180) * t113 + t130 * t116, t21 * qJD(3) - t64 * t74 + (t200 * qJDD(3) - t71 * t186) * pkin(3) - t137 + t3, t23 * qJD(3) + t64 * t71 + g(1) * t43 + g(2) * t41 + g(3) * t62 + (-qJDD(3) * t107 - t74 * t186) * pkin(3) - t4, (t19 - t21) * t74 + (-t18 + t23) * t71 + (-t107 * t38 - t200 * t39) * pkin(3), t3 * t172 - t64 * t177 + t18 * t21 - t19 * t23 + t4 * t215 + (-g(1) * (t108 * t192 + t113 * t67) - g(2) * (-t69 * t113 - t116 * t161) - t216) * pkin(3), t153 * t57 + t202, (t16 - t212) * t115 + (-t17 - t210) * t112, t153 * t65 + t203 - t209, t143 + t211, -t65 * t74, t127 * t112 - t115 * t218 + t145 * t74 + t96 * t17 - t21 * t55, t112 * t218 + t127 * t115 + t96 * t16 - t21 * t57 + t6 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t74 * qJD(3) + t146, (t90 - t71) * qJD(3) + t125, -t71 ^ 2 - t74 ^ 2, -t117 * t154 + t18 * t74 + t19 * t71 + t128 - t149, 0, 0, 0, 0, 0, t143 - t211, -t65 ^ 2 * t115 - t203 - t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t55, -t55 ^ 2 + t57 ^ 2, t16 + t212, -t17 + t210, t34, -t112 * t2 + t7 - t14 * t57 - g(1) * (-t112 * t43 + t115 * t70) - g(2) * (-t112 * t41 + t115 * t68) - g(3) * (-t112 * t62 - t173) + t221 * t6, -t115 * t2 - t112 * t8 + t14 * t55 - g(1) * (-t112 * t70 - t115 * t43) - g(2) * (-t112 * t68 - t115 * t41) - g(3) * (-t115 * t62 + t174) - t221 * t145;];
tau_reg = t5;
