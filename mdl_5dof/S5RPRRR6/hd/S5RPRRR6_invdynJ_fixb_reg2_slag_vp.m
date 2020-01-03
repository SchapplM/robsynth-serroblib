% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR6
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
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:46
% EndTime: 2019-12-31 19:01:51
% DurationCPUTime: 2.88s
% Computational Cost: add. (5233->402), mult. (11083->524), div. (0->0), fcn. (7471->14), ass. (0->220)
t144 = qJ(1) + pkin(9);
t136 = cos(t144);
t147 = qJ(3) + qJ(4);
t140 = sin(t147);
t234 = t136 * t140;
t135 = sin(t144);
t236 = t135 * t140;
t279 = g(1) * t234 + g(2) * t236;
t148 = sin(pkin(9));
t125 = pkin(1) * t148 + pkin(6);
t111 = t125 * qJDD(1);
t278 = qJD(2) * qJD(3) + t111;
t150 = sin(qJ(5));
t154 = cos(qJ(5));
t217 = qJD(3) + qJD(4);
t151 = sin(qJ(4));
t269 = cos(qJ(4));
t155 = cos(qJ(3));
t113 = t125 * qJD(1);
t200 = pkin(7) * qJD(1) + t113;
t152 = sin(qJ(3));
t226 = qJD(2) * t152;
t78 = t200 * t155 + t226;
t210 = t269 * t78;
t139 = t155 * qJD(2);
t77 = -t200 * t152 + t139;
t73 = qJD(3) * pkin(3) + t77;
t46 = t151 * t73 + t210;
t41 = t217 * pkin(8) + t46;
t205 = t269 * t155;
t227 = qJD(1) * t152;
t90 = -qJD(1) * t205 + t151 * t227;
t100 = t151 * t155 + t269 * t152;
t92 = t100 * qJD(1);
t149 = cos(pkin(9));
t126 = -pkin(1) * t149 - pkin(2);
t142 = t155 * pkin(3);
t275 = t126 - t142;
t93 = t275 * qJD(1);
t58 = pkin(4) * t90 - pkin(8) * t92 + t93;
t20 = -t150 * t41 + t154 * t58;
t21 = t150 * t58 + t154 * t41;
t183 = t150 * t21 + t154 * t20;
t222 = qJD(5) * t154;
t247 = t151 * t78;
t45 = t269 * t73 - t247;
t40 = -t217 * pkin(4) - t45;
t138 = t155 * qJDD(2);
t49 = qJDD(3) * pkin(3) + t138 + (-pkin(7) * qJDD(1) - t111) * t152 - t78 * qJD(3);
t221 = qJD(1) * qJD(3);
t203 = t152 * t221;
t218 = t155 * qJDD(1);
t211 = t152 * qJDD(2) + t278 * t155;
t225 = qJD(3) * t152;
t62 = -t113 * t225 + t211;
t54 = (-t203 + t218) * pkin(7) + t62;
t202 = t151 * t54 - t269 * t49;
t11 = -t46 * qJD(4) - t202;
t216 = qJDD(3) + qJDD(4);
t9 = -t216 * pkin(4) - t11;
t277 = t9 * t150 + t40 * t222;
t240 = pkin(1) * qJDD(1);
t189 = g(1) * t135 - g(2) * t136;
t276 = t189 * t140;
t141 = cos(t147);
t229 = t141 * pkin(4) + t140 * pkin(8);
t274 = -t150 * t20 + t154 * t21;
t273 = t100 * qJDD(1);
t173 = -t151 * t152 + t205;
t271 = t217 * qJD(1);
t56 = -t173 * t271 - t273;
t76 = t150 * t217 + t154 * t92;
t29 = qJD(5) * t76 - t150 * t56 - t154 * t216;
t89 = qJD(5) + t90;
t198 = t154 * t89;
t223 = qJD(5) * t150;
t219 = t152 * qJDD(1);
t184 = -qJDD(1) * t205 + t151 * t219;
t70 = t217 * t100;
t57 = qJD(1) * t70 + t184;
t53 = qJDD(5) + t57;
t243 = t154 * t53;
t197 = t151 * t217;
t204 = qJD(4) * t269;
t69 = -qJD(3) * t205 + t152 * t197 - t155 * t204;
t272 = -t100 * (t89 * t223 - t243) - t69 * t198;
t270 = t100 * t216 - t69 * t217;
t153 = sin(qJ(1));
t268 = pkin(1) * t153;
t130 = g(3) * t140;
t267 = g(3) * t141;
t266 = g(3) * t155;
t79 = pkin(3) * t203 + qJDD(1) * t275;
t19 = pkin(4) * t57 + pkin(8) * t56 + t79;
t224 = qJD(4) * t151;
t10 = t151 * t49 + t73 * t204 - t78 * t224 + t269 * t54;
t8 = t216 * pkin(8) + t10;
t2 = qJD(5) * t20 + t150 * t19 + t154 * t8;
t1 = t2 * t154;
t265 = t40 * t90;
t196 = t154 * t217;
t74 = t150 * t92 - t196;
t264 = t74 * t89;
t263 = t76 * t74;
t262 = t76 * t89;
t261 = t89 * t92;
t260 = t92 * t90;
t259 = pkin(7) + t125;
t242 = t154 * t69;
t244 = t154 * t29;
t258 = -t100 * t244 + t74 * t242;
t50 = t151 * t77 + t210;
t257 = t46 - t50;
t28 = -qJD(5) * t196 - t150 * t216 + t154 * t56 + t92 * t223;
t256 = t173 * t28 + t76 * t70;
t255 = -t100 * t57 + t69 * t90;
t252 = t150 * t28;
t251 = t150 * t53;
t250 = t150 * t69;
t249 = t150 * t74;
t248 = t150 * t89;
t241 = t154 * t76;
t239 = qJD(5) * t89;
t238 = t113 * t152;
t237 = t113 * t155;
t235 = t135 * t141;
t233 = t136 * t141;
t232 = t141 * t150;
t231 = t141 * t154;
t134 = t142 + pkin(2);
t156 = cos(qJ(1));
t143 = t156 * pkin(1);
t230 = t136 * t134 + t143;
t145 = t152 ^ 2;
t146 = t155 ^ 2;
t228 = t145 - t146;
t114 = qJD(1) * t126;
t112 = qJDD(1) * t126;
t215 = pkin(8) * t239;
t214 = t76 * t250;
t212 = pkin(3) * t225;
t132 = pkin(3) * t151 + pkin(8);
t209 = t132 * t239;
t159 = qJD(1) ^ 2;
t208 = t152 * t159 * t155;
t37 = t40 * t223;
t207 = g(1) * t233 + g(2) * t235 + t130;
t206 = -t9 - t267;
t201 = qJD(3) * t259;
t195 = pkin(3) * t204;
t194 = t155 * t203;
t193 = pkin(3) * t224 - t50;
t67 = pkin(4) * t92 + pkin(8) * t90;
t192 = t155 * t201;
t191 = -pkin(8) * t53 + t265;
t190 = g(1) * t136 + g(2) * t135;
t188 = g(1) * t153 - g(2) * t156;
t187 = t173 * t29 - t70 * t74;
t186 = t173 * t56 + t70 * t92;
t157 = -pkin(7) - pkin(6);
t185 = -t136 * t157 - t268;
t64 = -pkin(4) * t173 - pkin(8) * t100 + t275;
t94 = t259 * t152;
t95 = t259 * t155;
t66 = -t151 * t94 + t269 * t95;
t30 = -t150 * t66 + t154 * t64;
t31 = t150 * t64 + t154 * t66;
t181 = g(3) * t232 + t21 * t92 + t277;
t180 = t279 * t154 - t20 * t92 + t37;
t86 = t226 + t237;
t178 = -t183 * t90 + t1 - t207;
t177 = -qJD(5) * t58 + t130 - t8;
t176 = -t151 * t95 - t269 * t94;
t174 = t190 * t140;
t172 = -qJD(1) * t114 + t190;
t171 = t93 * t90 - t10 + t207;
t170 = 0.2e1 * t114 * qJD(3) - qJDD(3) * t125;
t18 = t154 * t19;
t3 = -qJD(5) * t21 - t150 * t8 + t18;
t169 = -qJD(5) * t183 - t3 * t150;
t168 = -t252 + (t241 + t249) * qJD(5);
t167 = -t93 * t92 - t202 - t267 + t279;
t166 = -t132 * t53 - t89 * t195 + t265;
t158 = qJD(3) ^ 2;
t165 = -t125 * t158 - 0.2e1 * t112 + t189;
t164 = (-t89 * t222 - t251) * t100 + t69 * t248;
t163 = t169 + t1;
t63 = -qJD(3) * t86 - t111 * t152 + t138;
t85 = t139 - t238;
t162 = -t63 * t152 + t62 * t155 + (-t152 * t86 - t155 * t85) * qJD(3);
t161 = -g(1) * (-pkin(4) * t234 + pkin(8) * t233) - g(2) * (-pkin(4) * t236 + pkin(8) * t235);
t133 = -t269 * pkin(3) - pkin(4);
t110 = qJDD(3) * t155 - t152 * t158;
t109 = qJDD(3) * t152 + t155 * t158;
t88 = t152 * t201;
t83 = t135 * t150 + t136 * t231;
t82 = t135 * t154 - t136 * t232;
t81 = -t135 * t231 + t136 * t150;
t80 = t135 * t232 + t136 * t154;
t61 = pkin(3) * t227 + t67;
t59 = -t90 ^ 2 + t92 ^ 2;
t55 = t173 * t216 - t70 * t217;
t51 = t269 * t77 - t247;
t36 = pkin(4) * t70 + pkin(8) * t69 + t212;
t35 = -t100 * t271 + t92 * t217 - t184;
t34 = t273 + (qJD(1) * t173 + t90) * t217;
t33 = t66 * qJD(4) - t151 * t88 + t269 * t192;
t32 = t176 * qJD(4) - t151 * t192 - t269 * t88;
t27 = t150 * t67 + t154 * t45;
t26 = -t150 * t45 + t154 * t67;
t25 = t150 * t61 + t154 * t51;
t24 = -t150 * t51 + t154 * t61;
t15 = t89 * t198 - t76 * t92 + t251;
t14 = -t89 ^ 2 * t150 + t74 * t92 + t243;
t13 = t74 * t248 - t244;
t12 = t76 * t198 - t252;
t6 = -qJD(5) * t31 - t150 * t32 + t154 * t36;
t5 = qJD(5) * t30 + t150 * t36 + t154 * t32;
t4 = (-t28 - t264) * t154 + (-t29 - t262) * t150;
t7 = [0, 0, 0, 0, 0, qJDD(1), t188, g(1) * t156 + g(2) * t153, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t149 * t240 + t189, -0.2e1 * t148 * t240 + t190, 0, (t188 + (t148 ^ 2 + t149 ^ 2) * t240) * pkin(1), qJDD(1) * t145 + 0.2e1 * t194, 0.2e1 * t152 * t218 - 0.2e1 * t228 * t221, t109, qJDD(1) * t146 - 0.2e1 * t194, t110, 0, t152 * t170 + t155 * t165, -t152 * t165 + t155 * t170, (t145 + t146) * t111 + t162 - t190, t112 * t126 - g(1) * (-pkin(2) * t135 + pkin(6) * t136 - t268) - g(2) * (pkin(2) * t136 + pkin(6) * t135 + t143) + t162 * t125, -t100 * t56 - t69 * t92, -t186 + t255, t270, -t173 * t57 + t70 * t90, t55, 0, t189 * t141 - t173 * t79 + t176 * t216 + t90 * t212 - t33 * t217 + t275 * t57 + t93 * t70, t79 * t100 + t92 * t212 - t66 * t216 - t32 * t217 - t275 * t56 - t93 * t69 - t276, t10 * t173 - t100 * t11 + t176 * t56 - t32 * t90 + t33 * t92 + t45 * t69 - t46 * t70 - t57 * t66 - t190, t10 * t66 + t46 * t32 + t11 * t176 - t45 * t33 + t79 * t275 + t93 * t212 - g(1) * (-t134 * t135 + t185) - g(2) * (-t135 * t157 + t230), -t69 * t241 + (-t154 * t28 - t76 * t223) * t100, t214 + (t252 + (-t241 + t249) * qJD(5)) * t100 + t258, t256 + t272, -t69 * t249 + (t150 * t29 + t74 * t222) * t100, t164 + t187, -t173 * t53 + t70 * t89, -g(1) * t81 - g(2) * t83 + t277 * t100 - t173 * t3 - t176 * t29 + t20 * t70 - t40 * t250 + t30 * t53 + t33 * t74 + t6 * t89, -t40 * t242 - g(1) * t80 - g(2) * t82 + t2 * t173 - t21 * t70 + t28 * t176 - t31 * t53 + t33 * t76 - t5 * t89 + (t9 * t154 - t37) * t100, t28 * t30 - t29 * t31 - t5 * t74 - t6 * t76 + t183 * t69 + t276 + (-qJD(5) * t274 - t150 * t2 - t154 * t3) * t100, t2 * t31 + t21 * t5 + t3 * t30 + t20 * t6 - t9 * t176 + t40 * t33 - g(1) * t185 - g(2) * (pkin(4) * t233 + pkin(8) * t234 + t230) + (-g(1) * (-t134 - t229) + g(2) * t157) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t110, -t109, 0, t152 * t62 + t155 * t63 - g(3) + (-t152 * t85 + t155 * t86) * qJD(3), 0, 0, 0, 0, 0, 0, t55, -t270, t186 + t255, t10 * t100 + t11 * t173 - t45 * t70 - t46 * t69 - g(3), 0, 0, 0, 0, 0, 0, t164 - t187, t256 - t272, t100 * t168 - t214 + t258, t100 * t163 - t173 * t9 - t274 * t69 + t40 * t70 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, t228 * t159, t219, t208, t218, qJDD(3), -t266 + t138 + (t86 - t237) * qJD(3) + (t172 - t278) * t152, g(3) * t152 + (t85 + t238) * qJD(3) + t172 * t155 - t211, 0, 0, t260, t59, t34, -t260, t35, t216, t50 * qJD(3) - t257 * qJD(4) + (-qJD(4) * t197 + t269 * t216 - t90 * t227) * pkin(3) + t167, t51 * t217 + (-t151 * t216 - t217 * t204 - t92 * t227) * pkin(3) + t171, t257 * t92 + (-t45 + t51) * t90 + (t269 * t56 - t151 * t57 + (t151 * t92 - t269 * t90) * qJD(4)) * pkin(3), t45 * t50 - t46 * t51 + (t269 * t11 - t266 + t10 * t151 + (-t151 * t45 + t269 * t46) * qJD(4) + (-qJD(1) * t93 + t190) * t152) * pkin(3), t12, t4, t15, t13, t14, -t261, t133 * t29 - t24 * t89 + t193 * t74 + (t206 - t209) * t154 + t166 * t150 + t180, -t133 * t28 + t25 * t89 + t193 * t76 + t166 * t154 + (-t174 + t209) * t150 + t181, t24 * t76 + t25 * t74 + (-t74 * t195 - t132 * t29 + (t132 * t76 - t20) * qJD(5)) * t154 + (t76 * t195 - t132 * t28 - t3 + (t132 * t74 - t21) * qJD(5)) * t150 + t178, t9 * t133 - t21 * t25 - t20 * t24 - t40 * t50 - g(3) * (t142 + t229) + (t190 * t152 + (t151 * t40 + t274 * t269) * qJD(4)) * pkin(3) + t163 * t132 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t59, t34, -t260, t35, t216, t46 * qJD(3) + t167, t45 * t217 + t171, 0, 0, t12, t4, t15, t13, t14, -t261, -pkin(4) * t29 - t26 * t89 - t46 * t74 + t191 * t150 + (t206 - t215) * t154 + t180, pkin(4) * t28 + t27 * t89 - t46 * t76 + t191 * t154 + (-t174 + t215) * t150 + t181, t26 * t76 + t27 * t74 + (t168 - t244) * pkin(8) + t169 + t178, -t9 * pkin(4) + t163 * pkin(8) - g(3) * t229 - t20 * t26 - t21 * t27 - t40 * t46 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, -t74 ^ 2 + t76 ^ 2, -t28 + t264, -t263, t262 - t29, t53, -g(1) * t82 + g(2) * t80 + t150 * t177 + t21 * t89 - t41 * t222 - t40 * t76 + t18, g(1) * t83 - g(2) * t81 + t20 * t89 + t40 * t74 + (qJD(5) * t41 - t19) * t150 + t177 * t154, 0, 0;];
tau_reg = t7;
