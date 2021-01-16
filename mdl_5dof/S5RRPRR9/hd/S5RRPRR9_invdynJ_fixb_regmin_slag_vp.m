% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:47:33
% EndTime: 2021-01-15 21:47:50
% DurationCPUTime: 3.63s
% Computational Cost: add. (3782->387), mult. (8934->529), div. (0->0), fcn. (6781->14), ass. (0->199)
t169 = cos(qJ(2));
t238 = cos(pkin(9));
t201 = t238 * t169;
t142 = qJD(1) * t201;
t161 = sin(pkin(9));
t165 = sin(qJ(2));
t223 = qJD(1) * t165;
t116 = t161 * t223 - t142;
t260 = qJD(4) + qJD(5);
t271 = t116 + t260;
t108 = qJD(4) + t116;
t145 = t161 * pkin(2) + pkin(7);
t157 = qJ(2) + pkin(9);
t150 = sin(t157);
t151 = cos(t157);
t166 = sin(qJ(1));
t170 = cos(qJ(1));
t194 = g(1) * t170 + g(2) * t166;
t179 = -g(3) * t151 + t194 * t150;
t162 = -qJ(3) - pkin(6);
t205 = qJD(2) * t162;
t180 = -t165 * qJD(3) + t169 * t205;
t209 = t162 * t165;
t76 = qJDD(2) * pkin(2) + t180 * qJD(1) + qJDD(1) * t209;
t115 = t169 * qJD(3) + t165 * t205;
t139 = t162 * t169;
t84 = t115 * qJD(1) - qJDD(1) * t139;
t30 = -t161 * t84 + t238 * t76;
t28 = -qJDD(2) * pkin(3) - t30;
t270 = -qJD(4) * t145 * t108 + t179 - t28;
t163 = sin(qJ(5));
t164 = sin(qJ(4));
t167 = cos(qJ(5));
t168 = cos(qJ(4));
t132 = t163 * t168 + t167 * t164;
t250 = t271 * t132;
t202 = t238 * t165;
t129 = t161 * t169 + t202;
t119 = t129 * qJD(1);
t218 = t168 * qJD(2);
t95 = t164 * t119 - t218;
t97 = t164 * qJD(2) + t168 * t119;
t189 = t163 * t95 - t167 * t97;
t39 = t163 * t97 + t167 * t95;
t269 = t189 * t39;
t222 = qJD(4) * t164;
t235 = t164 * t116;
t268 = t222 + t235;
t267 = t189 ^ 2 - t39 ^ 2;
t105 = qJD(5) + t108;
t219 = qJD(5) * t167;
t220 = qJD(5) * t163;
t217 = qJD(1) * qJD(2);
t208 = t165 * t217;
t175 = t129 * qJDD(1) - t161 * t208;
t83 = qJD(2) * t142 + t175;
t35 = qJD(4) * t218 + t164 * qJDD(2) - t119 * t222 + t168 * t83;
t36 = t97 * qJD(4) - t168 * qJDD(2) + t164 * t83;
t7 = -t163 * t36 + t167 * t35 - t95 * t219 - t97 * t220;
t266 = t39 * t105 + t7;
t160 = qJ(4) + qJ(5);
t156 = cos(t160);
t227 = t170 * t156;
t155 = sin(t160);
t233 = t166 * t155;
t101 = t151 * t227 + t233;
t149 = t169 * pkin(2) + pkin(1);
t137 = -t149 * qJD(1) + qJD(3);
t48 = t116 * pkin(3) - t119 * pkin(7) + t137;
t133 = qJD(1) * t209;
t249 = qJD(2) * pkin(2);
t125 = t133 + t249;
t134 = qJD(1) * t139;
t203 = t238 * t134;
t81 = t161 * t125 - t203;
t70 = qJD(2) * pkin(7) + t81;
t23 = t164 * t48 + t168 * t70;
t17 = -t95 * pkin(8) + t23;
t14 = t17 * t220;
t257 = g(3) * t150;
t122 = t161 * t134;
t80 = t238 * t125 + t122;
t69 = -qJD(2) * pkin(3) - t80;
t37 = t95 * pkin(4) + t69;
t228 = t170 * t155;
t232 = t166 * t156;
t99 = -t151 * t232 + t228;
t265 = g(1) * t101 - g(2) * t99 + t156 * t257 + t37 * t39 + t14;
t100 = -t151 * t228 + t232;
t107 = pkin(2) * t208 - t149 * qJDD(1) + qJDD(3);
t118 = t129 * qJD(2);
t216 = t165 * qJDD(1);
t190 = -qJDD(1) * t201 + t161 * t216;
t82 = qJD(1) * t118 + t190;
t24 = t82 * pkin(3) - t83 * pkin(7) + t107;
t21 = t168 * t24;
t31 = t161 * t76 + t238 * t84;
t29 = qJDD(2) * pkin(7) + t31;
t77 = qJDD(4) + t82;
t2 = t77 * pkin(4) - t35 * pkin(8) - t23 * qJD(4) - t164 * t29 + t21;
t221 = qJD(4) * t168;
t186 = t164 * t24 + t168 * t29 + t48 * t221 - t70 * t222;
t3 = -t36 * pkin(8) + t186;
t213 = -t163 * t3 + t167 * t2;
t22 = -t164 * t70 + t168 * t48;
t16 = -t97 * pkin(8) + t22;
t11 = t108 * pkin(4) + t16;
t244 = t167 * t17;
t5 = t163 * t11 + t244;
t98 = t151 * t233 + t227;
t264 = -g(1) * t100 + g(2) * t98 - t5 * qJD(5) + t155 * t257 + t37 * t189 + t213;
t174 = t189 * qJD(5) - t163 * t35 - t167 * t36;
t263 = -t105 * t189 + t174;
t65 = t132 * t129;
t181 = -t161 * t165 + t201;
t121 = t181 * qJD(2);
t229 = t168 * t121;
t261 = -t129 * t222 + t229;
t131 = t163 * t164 - t167 * t168;
t251 = t271 * t131;
t72 = qJDD(5) + t77;
t259 = t251 * t105 - t132 * t72;
t258 = pkin(2) * t165;
t255 = g(3) * t169;
t254 = pkin(8) + t145;
t59 = pkin(2) * t223 + t119 * pkin(3) + t116 * pkin(7);
t86 = t238 * t133 + t122;
t253 = t164 * t59 + t168 * t86;
t79 = -pkin(3) * t181 - t129 * pkin(7) - t149;
t94 = -t238 * t139 + t161 * t209;
t87 = t168 * t94;
t252 = t164 * t79 + t87;
t248 = t119 * t39;
t247 = t119 * t95;
t245 = t164 * t77;
t243 = t35 * t164;
t242 = t189 * t119;
t241 = t95 * t108;
t240 = t97 * t108;
t239 = t97 * t119;
t237 = t129 * t164;
t236 = t129 * t168;
t234 = t164 * t121;
t231 = t166 * t164;
t230 = t166 * t168;
t226 = t170 * t164;
t225 = t170 * t168;
t158 = t165 ^ 2;
t224 = -t169 ^ 2 + t158;
t215 = t169 * qJDD(1);
t214 = t165 * t249;
t210 = t129 * t221;
t207 = qJD(5) * t11 + t3;
t204 = qJD(4) * t254;
t200 = -qJD(4) * t48 - t29;
t57 = t161 * t115 - t238 * t180;
t85 = t161 * t133 - t203;
t93 = -t161 * t139 - t162 * t202;
t198 = t108 * t168;
t197 = -t250 * t105 - t131 * t72;
t196 = t268 * pkin(4) - t85;
t146 = -t238 * pkin(2) - pkin(3);
t195 = -t70 * t221 + t21;
t193 = g(1) * t166 - g(2) * t170;
t126 = t254 * t164;
t192 = pkin(8) * t235 + qJD(5) * t126 + t164 * t204 + t253;
t127 = t254 * t168;
t52 = t168 * t59;
t191 = t119 * pkin(4) + qJD(5) * t127 - t164 * t86 + t52 + (pkin(8) * t116 + t204) * t168;
t188 = -t268 * t108 + t168 * t77;
t187 = -0.2e1 * pkin(1) * t217 - pkin(6) * qJDD(2);
t58 = t238 * t115 + t161 * t180;
t60 = t118 * pkin(3) - t121 * pkin(7) + t214;
t185 = t164 * t60 + t168 * t58 + t79 * t221 - t94 * t222;
t184 = t210 + t234;
t182 = t108 * t69 - t145 * t77;
t171 = qJD(2) ^ 2;
t177 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t171 + t193;
t172 = qJD(1) ^ 2;
t176 = pkin(1) * t172 - pkin(6) * qJDD(1) + t194;
t138 = -t168 * pkin(4) + t146;
t112 = t151 * t225 + t231;
t111 = -t151 * t226 + t230;
t110 = -t151 * t230 + t226;
t109 = t151 * t231 + t225;
t68 = t168 * t79;
t66 = t131 * t129;
t56 = pkin(4) * t237 + t93;
t53 = t168 * t60;
t27 = t184 * pkin(4) + t57;
t25 = -pkin(8) * t237 + t252;
t18 = -pkin(4) * t181 - pkin(8) * t236 - t164 * t94 + t68;
t13 = -t220 * t237 + (t260 * t236 + t234) * t167 + t261 * t163;
t12 = -t131 * t121 - t260 * t65;
t10 = t36 * pkin(4) + t28;
t9 = -t184 * pkin(8) + t185;
t6 = -pkin(8) * t229 + t118 * pkin(4) - t164 * t58 + t53 + (-t87 + (pkin(8) * t129 - t79) * t164) * qJD(4);
t4 = t167 * t11 - t163 * t17;
t1 = [qJDD(1), t193, t194, t158 * qJDD(1) + 0.2e1 * t169 * t208, 0.2e1 * t165 * t215 - 0.2e1 * t224 * t217, qJDD(2) * t165 + t171 * t169, qJDD(2) * t169 - t171 * t165, 0, t187 * t165 + t177 * t169, -t177 * t165 + t187 * t169, -t93 * qJDD(2) - t107 * t181 + t137 * t118 - t149 * t82 + t193 * t151 + (t116 * t258 - t57) * qJD(2), -t94 * qJDD(2) + t107 * t129 + t137 * t121 - t149 * t83 - t193 * t150 + (t119 * t258 - t58) * qJD(2), -t58 * t116 - t81 * t118 + t57 * t119 - t80 * t121 - t30 * t129 + t181 * t31 - t94 * t82 + t93 * t83 - t194, t31 * t94 + t81 * t58 - t30 * t93 - t80 * t57 - t107 * t149 + t137 * t214 - g(1) * (-t166 * t149 - t162 * t170) - g(2) * (t170 * t149 - t166 * t162), t97 * t229 + (t35 * t168 - t97 * t222) * t129, (-t164 * t97 - t168 * t95) * t121 + (-t243 - t168 * t36 + (t164 * t95 - t168 * t97) * qJD(4)) * t129, t261 * t108 + t97 * t118 - t181 * t35 + t77 * t236, -t108 * t184 - t95 * t118 + t181 * t36 - t77 * t237, t108 * t118 - t181 * t77, (-t94 * t221 + t53) * t108 + t68 * t77 - t195 * t181 + t22 * t118 + t57 * t95 + t93 * t36 + t69 * t210 - g(1) * t110 - g(2) * t112 + ((-qJD(4) * t79 - t58) * t108 - t94 * t77 - t200 * t181 + t28 * t129 + t69 * t121) * t164, -t185 * t108 - t252 * t77 + t186 * t181 - t23 * t118 + t57 * t97 + t93 * t35 + t69 * t229 - g(1) * t109 - g(2) * t111 + (t28 * t168 - t69 * t222) * t129, -t12 * t189 - t7 * t66, -t12 * t39 + t13 * t189 - t174 * t66 - t7 * t65, t12 * t105 - t118 * t189 - t181 * t7 - t66 * t72, -t13 * t105 - t39 * t118 - t174 * t181 - t65 * t72, t105 * t118 - t181 * t72, (-t163 * t9 + t167 * t6) * t105 + (-t163 * t25 + t167 * t18) * t72 - t213 * t181 + t4 * t118 + t27 * t39 - t56 * t174 + t10 * t65 + t37 * t13 - g(1) * t99 - g(2) * t101 + ((-t163 * t18 - t167 * t25) * t105 + t5 * t181) * qJD(5), -g(1) * t98 - g(2) * t100 - t10 * t66 - t5 * t118 + t37 * t12 - t14 * t181 - t27 * t189 + t56 * t7 + (-(-qJD(5) * t25 + t6) * t105 - t18 * t72 + t2 * t181) * t163 + (-(qJD(5) * t18 + t9) * t105 - t25 * t72 + t207 * t181) * t167; 0, 0, 0, -t165 * t172 * t169, t224 * t172, t216, t215, qJDD(2), t165 * t176 - t255, g(3) * t165 + t169 * t176, t85 * qJD(2) - t137 * t119 + (t238 * qJDD(2) - t116 * t223) * pkin(2) + t179 + t30, t257 + t86 * qJD(2) + t137 * t116 + t194 * t151 + (-qJDD(2) * t161 - t119 * t223) * pkin(2) - t31, (t81 - t85) * t119 + (-t80 + t86) * t116 + (-t161 * t82 - t238 * t83) * pkin(2), t80 * t85 - t81 * t86 + (t238 * t30 - t255 + t161 * t31 + (-qJD(1) * t137 + t194) * t165) * pkin(2), t198 * t97 + t243, (t35 - t241) * t168 + (-t36 - t240) * t164, t108 * t198 - t239 + t245, t188 + t247, -t108 * t119, -t52 * t108 - t22 * t119 + t146 * t36 - t85 * t95 + (t86 * t108 + t182) * t164 + t270 * t168, t253 * t108 + t23 * t119 + t146 * t35 - t270 * t164 + t182 * t168 - t85 * t97, t7 * t132 + t189 * t251, -t7 * t131 + t132 * t174 + t189 * t250 + t251 * t39, t242 - t259, t197 + t248, -t105 * t119, (-t167 * t126 - t163 * t127) * t72 - t138 * t174 + t10 * t131 - t4 * t119 + t196 * t39 + t250 * t37 + (t163 * t192 - t167 * t191) * t105 + t179 * t156, -(-t163 * t126 + t167 * t127) * t72 + t138 * t7 + t10 * t132 + t5 * t119 - t196 * t189 - t251 * t37 + (t163 * t191 + t167 * t192) * t105 - t179 * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t119 * qJD(2) + t190, (t142 - t116) * qJD(2) + t175, -t116 ^ 2 - t119 ^ 2, t81 * t116 + t80 * t119 + t107 - t193, 0, 0, 0, 0, 0, t188 - t247, -t108 ^ 2 * t168 - t239 - t245, 0, 0, 0, 0, 0, t197 - t248, t242 + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t95, -t95 ^ 2 + t97 ^ 2, t35 + t241, t240 - t36, t77, -g(1) * t111 + g(2) * t109 + t23 * t108 - t69 * t97 + (t200 + t257) * t164 + t195, g(1) * t112 - g(2) * t110 + t22 * t108 + t168 * t257 + t69 * t95 - t186, -t269, t267, t266, t263, t72, -(-t163 * t16 - t244) * t105 + (-t105 * t220 + t167 * t72 - t97 * t39) * pkin(4) + t264, (-t17 * t105 - t2) * t163 + (t16 * t105 - t207) * t167 + (-t105 * t219 - t163 * t72 + t189 * t97) * pkin(4) + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t269, t267, t266, t263, t72, t5 * t105 + t264, t4 * t105 - t163 * t2 - t167 * t207 + t265;];
tau_reg = t1;
