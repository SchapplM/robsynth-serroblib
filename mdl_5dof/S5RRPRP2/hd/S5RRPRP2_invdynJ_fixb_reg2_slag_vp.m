% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:47
% EndTime: 2019-12-31 19:49:50
% DurationCPUTime: 1.27s
% Computational Cost: add. (1902->253), mult. (3120->281), div. (0->0), fcn. (1762->12), ass. (0->159)
t115 = sin(qJ(4));
t110 = t115 ^ 2;
t118 = cos(qJ(4));
t111 = t118 ^ 2;
t222 = t110 + t111;
t108 = qJDD(1) + qJDD(2);
t185 = t111 * t108;
t186 = t110 * t108;
t221 = t185 + t186;
t113 = sin(pkin(8));
t114 = cos(pkin(8));
t119 = cos(qJ(2));
t198 = pkin(1) * qJD(2);
t163 = qJD(1) * t198;
t116 = sin(qJ(2));
t171 = qJDD(1) * t116;
t218 = pkin(1) * t171 + t119 * t163;
t204 = t119 * pkin(1);
t98 = qJDD(1) * t204;
t43 = t108 * pkin(2) - t116 * t163 + t98;
t19 = t113 * t43 + t218 * t114;
t220 = t108 * pkin(7) + qJD(3) * qJD(4) + t19;
t109 = qJD(1) + qJD(2);
t199 = pkin(1) * qJD(1);
t165 = t116 * t199;
t164 = t119 * t199;
t61 = t109 * pkin(2) + t164;
t39 = t113 * t61 + t114 * t165;
t30 = t109 * pkin(7) + t39;
t197 = t115 * t30;
t21 = t118 * qJD(3) - t197;
t219 = qJD(5) - t21;
t143 = t118 * pkin(4) + t115 * qJ(5);
t166 = t115 * qJDD(3) + t220 * t118;
t170 = qJDD(4) * qJ(5);
t3 = t170 + (qJD(5) - t197) * qJD(4) + t166;
t176 = qJD(4) * t118;
t156 = -t118 * qJDD(3) + t220 * t115 + t30 * t176;
t188 = qJDD(4) * pkin(4);
t215 = qJDD(5) - t188;
t4 = t156 + t215;
t217 = t4 * t115 + t3 * t118;
t22 = t115 * qJD(3) + t118 * t30;
t177 = qJD(4) * t115;
t6 = -t30 * t177 + t166;
t124 = t6 * t118 + (-t115 * t22 - t118 * t21) * qJD(4) + t156 * t115;
t112 = qJ(1) + qJ(2);
t102 = pkin(8) + t112;
t87 = sin(t102);
t88 = cos(t102);
t201 = g(1) * t88 + g(2) * t87;
t17 = -qJD(4) * pkin(4) + t219;
t20 = qJD(4) * qJ(5) + t22;
t104 = sin(t112);
t105 = cos(t112);
t216 = g(1) * t104 - g(2) * t105;
t77 = t113 * t165;
t50 = t114 * t164 - t77;
t190 = t50 * t109;
t207 = t113 * pkin(2);
t89 = pkin(7) + t207;
t214 = -t222 * t190 + t221 * t89;
t184 = t113 * t116;
t51 = (t114 * t119 - t184) * t198;
t189 = t51 * t109;
t183 = t114 * t116;
t97 = pkin(2) + t204;
t53 = pkin(1) * t183 + t113 * t97;
t47 = pkin(7) + t53;
t213 = t222 * t189 + t221 * t47;
t178 = qJD(4) * t109;
t180 = -t110 + t111;
t181 = t118 * t108;
t212 = 0.2e1 * t115 * t181 + 0.2e1 * t180 * t178;
t82 = g(1) * t87;
t211 = g(2) * t88;
t210 = t87 * pkin(3);
t209 = pkin(2) * t104;
t206 = t114 * pkin(2);
t134 = pkin(1) * (t113 * t119 + t183);
t48 = qJD(1) * t134;
t173 = t115 * qJD(5);
t54 = pkin(4) * t177 - qJ(5) * t176 - t173;
t203 = t54 - t48;
t195 = t115 * t88;
t196 = t115 * t87;
t202 = g(1) * t196 - g(2) * t195;
t200 = g(1) * t105 + g(2) * t104;
t121 = qJD(4) ^ 2;
t194 = t121 * t47;
t193 = t121 * t89;
t192 = t48 * t109;
t49 = qJD(2) * t134;
t191 = t49 * t109;
t187 = t109 * t115;
t175 = qJDD(4) * t47;
t174 = qJDD(4) * t89;
t69 = t118 * t82;
t169 = t118 * t192 + t50 * t177 + t69;
t168 = t218 * t113 - t114 * t43;
t96 = pkin(2) * t105;
t167 = t88 * pkin(3) + t87 * pkin(7) + t96;
t13 = -t108 * pkin(3) + t168;
t161 = -t13 - t211;
t79 = t88 * pkin(7);
t160 = t79 - t209;
t137 = -pkin(3) - t143;
t52 = -pkin(1) * t184 + t114 * t97;
t31 = t137 - t52;
t159 = t109 * t31 - t51;
t38 = t114 * t61 - t77;
t158 = t21 + t197;
t29 = -t109 * pkin(3) - t38;
t157 = t13 * t115 + t29 * t176 - t202;
t154 = qJD(1) * (-qJD(2) + t109);
t153 = qJD(2) * (-qJD(1) - t109);
t151 = t98 + t216;
t150 = t176 * t187;
t149 = t143 * t88 + t167;
t148 = -t19 + t201;
t117 = sin(qJ(1));
t145 = -t117 * pkin(1) - t209;
t120 = cos(qJ(1));
t144 = g(1) * t117 - g(2) * t120;
t142 = pkin(4) * t115 - qJ(5) * t118;
t90 = -pkin(3) - t206;
t141 = t108 * t90 + t193;
t140 = t17 * t115 + t20 * t118;
t139 = t21 * t115 - t22 * t118;
t138 = t145 + t79;
t57 = t137 - t206;
t8 = (t142 * qJD(4) - t173) * t109 + t137 * t108 + t168;
t136 = -t108 * t57 - t193 - t8;
t135 = -t168 + t82 - t211;
t46 = -pkin(3) - t52;
t133 = t108 * t46 + t191 + t194;
t132 = g(1) * t195 + g(2) * t196 - g(3) * t118 - t156;
t131 = t17 * t176 - t20 * t177 - t201 + t217;
t130 = t137 * t82;
t23 = t49 + t54;
t129 = -t108 * t31 - t109 * t23 - t194 - t8;
t128 = -t175 + (t109 * t46 - t51) * qJD(4);
t126 = t22 * qJD(4) + t132;
t125 = (-t115 * t20 + t118 * t17) * qJD(4) + t217;
t123 = -t201 + t124;
t107 = t109 ^ 2;
t106 = t120 * pkin(1);
t86 = t115 * t108;
t73 = t115 * t107 * t118;
t64 = qJDD(4) * t118 - t121 * t115;
t63 = qJDD(4) * t115 + t121 * t118;
t58 = t180 * t107;
t55 = t142 * t109;
t45 = -0.2e1 * t150 + t185;
t44 = 0.2e1 * t150 + t186;
t24 = t29 * t177;
t16 = t137 * t109 - t38;
t12 = t16 * t177;
t1 = [0, 0, 0, 0, 0, qJDD(1), t144, g(1) * t120 + g(2) * t117, 0, 0, 0, 0, 0, 0, 0, t108, (t108 * t119 + t116 * t153) * pkin(1) + t151, ((-qJDD(1) - t108) * t116 + t119 * t153) * pkin(1) + t200, 0, (t144 + (t116 ^ 2 + t119 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t108, t52 * t108 + t135 - t191, -t53 * t108 + t148 - t189, 0, t19 * t53 + t39 * t51 - t168 * t52 - t38 * t49 - g(1) * t145 - g(2) * (t96 + t106), t44, t212, t63, t45, t64, 0, t24 + t69 + t128 * t115 + (-t133 + t161) * t118, t133 * t115 + t128 * t118 + t157, t123 + t213, t13 * t46 + t29 * t49 - g(1) * (t138 - t210) - g(2) * (t106 + t167) - t139 * t51 + t124 * t47, t44, t63, -t212, 0, -t64, t45, t12 + t69 + (qJD(4) * t159 - t175) * t115 + (t129 - t211) * t118, t131 + t213, (t175 + (-t159 - t16) * qJD(4)) * t118 + t129 * t115 + t202, t8 * t31 + t16 * t23 - g(1) * t138 - g(2) * (t106 + t149) - t130 + t140 * t51 + t125 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t116 * pkin(1) * t154 + t151, (t119 * t154 - t171) * pkin(1) + t200, 0, 0, 0, 0, 0, 0, 0, t108, t108 * t206 + t135 + t192, -t108 * t207 + t148 + t190, 0, t38 * t48 - t39 * t50 + (t113 * t19 - t114 * t168 + t216) * pkin(2), t44, t212, t63, t45, t64, 0, t24 + (t90 * t178 - t174) * t115 + (-t141 + t161) * t118 + t169, (-t174 + (t109 * t90 + t50) * qJD(4)) * t118 + (t141 - t192) * t115 + t157, t123 + t214, t13 * t90 - t29 * t48 - g(1) * (t160 - t210) - g(2) * t167 + t139 * t50 + t124 * t89, t44, t63, -t212, 0, -t64, t45, t12 + (t178 * t57 - t174) * t115 + (-t109 * t54 + t136 - t211) * t118 + t169, t131 + t214, (t174 + (-t109 * t57 - t16 - t50) * qJD(4)) * t118 + (-t203 * t109 + t136) * t115 + t202, -g(1) * t160 - g(2) * t149 + t125 * t89 - t140 * t50 + t203 * t16 + t8 * t57 - t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, 0, 0, 0, 0, t64, -t63, 0, -t139 * qJD(4) + t6 * t115 - t118 * t156 - g(3), 0, 0, 0, 0, 0, 0, t64, 0, t63, qJD(4) * t140 + t3 * t115 - t4 * t118 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t58, t86, t73, t181, qJDD(4), -t29 * t187 + t126, g(3) * t115 + t158 * qJD(4) + (-t109 * t29 + t201) * t118 - t166, 0, 0, -t73, t86, t58, qJDD(4), -t181, t73, 0.2e1 * t188 - qJDD(5) + (-t115 * t16 + t118 * t55) * t109 + t126, -t142 * t108, 0.2e1 * t170 + (t109 * t55 - g(3)) * t115 + (t109 * t16 - t201) * t118 + (0.2e1 * qJD(5) - t158) * qJD(4) + t166, -t4 * pkin(4) - g(3) * t143 + t3 * qJ(5) + t201 * t142 - t16 * t55 - t17 * t22 + t219 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t73, t86, -t110 * t107 - t121, -t20 * qJD(4) + t16 * t187 - t132 + t215;];
tau_reg = t1;
