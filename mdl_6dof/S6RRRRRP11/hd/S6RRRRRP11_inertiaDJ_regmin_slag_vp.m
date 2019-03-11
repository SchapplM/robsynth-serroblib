% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRRP11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:53:38
% EndTime: 2019-03-10 02:54:00
% DurationCPUTime: 8.02s
% Computational Cost: add. (12990->549), mult. (39223->1000), div. (0->0), fcn. (39588->12), ass. (0->230)
t184 = cos(qJ(2));
t264 = qJD(2) * t184;
t176 = cos(pkin(6));
t291 = pkin(1) * t176;
t157 = t264 * t291;
t180 = sin(qJ(2));
t174 = sin(pkin(6));
t175 = cos(pkin(7));
t218 = t174 * (-pkin(10) * t175 - pkin(9));
t209 = t180 * t218;
t102 = qJD(2) * t209 + t157;
t179 = sin(qJ(3));
t101 = (pkin(1) * t184 + pkin(2)) * t176 + t209;
t173 = sin(pkin(7));
t290 = pkin(10) * t173;
t119 = (-pkin(2) * t184 - t180 * t290 - pkin(1)) * t174;
t211 = t101 * t175 + t119 * t173;
t183 = cos(qJ(3));
t249 = t180 * t291;
t280 = t174 * t184;
t196 = -pkin(9) * t280 - t249;
t94 = (t173 * t176 + t175 * t280) * pkin(10) - t196;
t91 = t183 * t94;
t298 = -(t179 * t211 + t91) * qJD(3) - t179 * t102;
t178 = sin(qJ(4));
t297 = -0.4e1 * t178;
t182 = cos(qJ(4));
t256 = qJD(4) * t183;
t236 = t182 * t256;
t263 = qJD(3) * t179;
t296 = (t178 * t263 - t236) * t173;
t216 = -pkin(4) * t182 - pkin(12) * t178;
t151 = -pkin(3) + t216;
t181 = cos(qJ(5));
t272 = t181 * t182;
t163 = pkin(11) * t272;
t177 = sin(qJ(5));
t125 = t151 * t177 + t163;
t271 = t183 * t184;
t275 = t179 * t180;
t295 = t175 * t271 - t275;
t171 = t181 ^ 2;
t267 = t177 ^ 2 - t171;
t226 = t267 * qJD(5);
t293 = -t179 * t94 + t183 * t211;
t172 = t182 ^ 2;
t292 = 0.2e1 * t173;
t289 = pkin(11) * t177;
t288 = pkin(12) * t182;
t287 = -qJ(6) - pkin(12);
t70 = -t101 * t173 + t119 * t175;
t281 = t173 * t183;
t98 = -t174 * t295 - t176 * t281;
t273 = t180 * t183;
t274 = t179 * t184;
t200 = t175 * t274 + t273;
t282 = t173 * t179;
t99 = t174 * t200 + t176 * t282;
t48 = pkin(3) * t98 - pkin(11) * t99 + t70;
t137 = t173 * t280 - t175 * t176;
t279 = t175 * t179;
t247 = t101 * t279 + t119 * t282 + t91;
t53 = -pkin(11) * t137 + t247;
t27 = t178 * t48 + t182 * t53;
t25 = pkin(12) * t98 + t27;
t52 = t137 * pkin(3) - t293;
t79 = t137 * t182 + t178 * t99;
t80 = -t137 * t178 + t182 * t99;
t31 = t79 * pkin(4) - t80 * pkin(12) + t52;
t13 = t177 * t31 + t181 * t25;
t159 = pkin(10) * t282;
t127 = t159 + (-pkin(2) * t183 - pkin(3)) * t175;
t138 = -t175 * t182 + t178 * t282;
t139 = t175 * t178 + t182 * t282;
t81 = t138 * pkin(4) - t139 * pkin(12) + t127;
t162 = pkin(2) * t279;
t268 = pkin(10) * t281 + t162;
t128 = pkin(11) * t175 + t268;
t217 = -pkin(3) * t183 - pkin(11) * t179;
t129 = (-pkin(2) + t217) * t173;
t87 = t128 * t182 + t129 * t178;
t83 = -pkin(12) * t281 + t87;
t55 = t177 * t81 + t181 * t83;
t255 = qJD(5) * t177;
t265 = qJD(2) * t180;
t231 = t174 * t265;
t221 = t173 * t231;
t261 = qJD(3) * t183;
t232 = t173 * t261;
t74 = t176 * t232 + (t295 * qJD(3) + (-t175 * t275 + t271) * qJD(2)) * t174;
t44 = -qJD(4) * t79 + t178 * t221 + t74 * t182;
t233 = t173 * t263;
t73 = t176 * t233 + (t200 * qJD(3) + (t175 * t273 + t274) * qJD(2)) * t174;
t22 = t73 * t177 - t80 * t255 + (qJD(5) * t98 + t44) * t181;
t286 = t22 * t177;
t285 = t22 * t181;
t105 = t139 * t177 + t181 * t281;
t190 = qJD(4) * t138;
t219 = t182 * t232;
t187 = -t190 + t219;
t66 = qJD(5) * t105 - t177 * t233 - t181 * t187;
t284 = t66 * t177;
t283 = qJ(6) * t178;
t278 = t175 * t183;
t277 = t178 * t181;
t215 = pkin(4) * t178 - t288;
t199 = t215 * qJD(4);
t254 = qJD(5) * t181;
t270 = -t151 * t254 - t177 * t199;
t259 = qJD(4) * t178;
t238 = t177 * t259;
t269 = pkin(11) * t238 + t181 * t199;
t170 = t178 ^ 2;
t266 = t170 - t172;
t262 = qJD(3) * t182;
t260 = qJD(4) * t177;
t258 = qJD(4) * t181;
t257 = qJD(4) * t182;
t253 = qJD(5) * t182;
t252 = t181 * qJD(6);
t251 = -0.2e1 * pkin(3) * qJD(4);
t250 = -0.2e1 * pkin(4) * qJD(5);
t248 = t182 * t289;
t246 = pkin(11) * t257;
t245 = pkin(5) * t255;
t243 = t177 * t281;
t168 = t174 ^ 2;
t242 = t168 * t264;
t167 = t173 ^ 2;
t241 = t167 * t261;
t240 = t175 * t261;
t239 = t175 * t257;
t237 = t181 * t257;
t235 = t177 * t253;
t234 = t181 * t253;
t230 = t174 * t264;
t229 = t177 * t254;
t228 = t178 * t257;
t12 = -t177 * t25 + t181 * t31;
t54 = -t177 * t83 + t181 * t81;
t26 = -t178 * t53 + t182 * t48;
t227 = qJD(5) * t287;
t103 = (t184 * t218 - t249) * qJD(2);
t123 = (pkin(2) * t180 - t184 * t290) * t174 * qJD(2);
t75 = -t103 * t173 + t123 * t175;
t86 = -t128 * t178 + t129 * t182;
t225 = t266 * qJD(4);
t224 = 0.2e1 * t228;
t223 = t167 * t231;
t222 = t179 * t241;
t220 = t177 * t237;
t82 = pkin(4) * t281 - t86;
t58 = t177 * t80 - t181 * t98;
t11 = -qJ(6) * t58 + t13;
t59 = t177 * t98 + t181 * t80;
t8 = pkin(5) * t79 - qJ(6) * t59 + t12;
t214 = -t11 * t177 - t181 * t8;
t106 = t139 * t181 - t243;
t41 = pkin(5) * t138 - qJ(6) * t106 + t54;
t42 = -qJ(6) * t105 + t55;
t213 = -t177 * t42 - t181 * t41;
t212 = -t177 * t59 - t181 * t58;
t36 = -t101 * t240 - t102 * t183 - t103 * t279 - t119 * t232 - t123 * t282 + t263 * t94;
t34 = pkin(11) * t221 - t36;
t39 = pkin(3) * t73 - pkin(11) * t74 + t75;
t10 = -t178 * t34 + t182 * t39 - t257 * t53 - t259 * t48;
t210 = -t105 * t181 - t106 * t177;
t24 = -pkin(4) * t98 - t26;
t130 = (pkin(3) * t179 - pkin(11) * t183) * t173 * qJD(3);
t131 = -pkin(2) * t240 + pkin(10) * t233;
t61 = -t128 * t257 - t129 * t259 + t130 * t182 + t131 * t178;
t7 = -pkin(4) * t73 - t10;
t208 = t7 * t177 + t24 * t254;
t207 = -t7 * t181 + t24 * t255;
t206 = t178 * t73 + t257 * t98;
t205 = -t182 * t73 + t259 * t98;
t43 = qJD(4) * t80 + t74 * t178 - t182 * t221;
t204 = t177 * t43 + t254 * t79;
t203 = -t181 * t43 + t255 * t79;
t57 = -pkin(4) * t233 - t61;
t202 = t57 * t177 + t254 * t82;
t201 = -t57 * t181 + t255 * t82;
t35 = -t103 * t278 + (-pkin(3) * t231 - t123 * t183) * t173 - t298;
t185 = t43 * pkin(4) - t44 * pkin(12) + t35;
t9 = -t178 * t39 - t182 * t34 - t257 * t48 + t259 * t53;
t189 = pkin(12) * t73 - t9;
t3 = -t177 * t185 - t181 * t189 + t25 * t255 - t254 * t31;
t104 = qJD(4) * t139 + t178 * t232;
t186 = t104 * pkin(4) + pkin(12) * t190 + (t162 + (pkin(10) - t288) * t281) * qJD(3);
t60 = t128 * t259 - t129 * t257 - t130 * t178 + t131 * t182;
t188 = pkin(12) * t233 - t60;
t19 = -t177 * t186 - t181 * t188 - t254 * t81 + t255 * t83;
t198 = t104 * t177 + t138 * t254;
t197 = -t104 * t181 + t138 * t255;
t195 = -t179 * t259 + t182 * t261;
t193 = -t178 * t255 + t237;
t192 = t178 * t258 + t235;
t191 = t177 * t257 + t178 * t254;
t4 = -qJD(5) * t13 - t177 * t189 + t181 * t185;
t20 = -qJD(5) * t55 - t177 * t188 + t181 * t186;
t166 = -pkin(5) * t181 - pkin(4);
t153 = t287 * t181;
t152 = t287 * t177;
t147 = (pkin(5) * t177 + pkin(11)) * t178;
t145 = t181 * t151;
t136 = -t177 * qJD(6) + t181 * t227;
t135 = t177 * t227 + t252;
t134 = t196 * qJD(2);
t133 = pkin(9) * t231 - t157;
t132 = t268 * qJD(3);
t124 = t145 - t248;
t118 = pkin(5) * t191 + t246;
t108 = -t177 * t283 + t125;
t97 = -qJ(6) * t277 + t145 + (-pkin(5) - t289) * t182;
t85 = -qJD(5) * t125 + t269;
t84 = pkin(11) * t192 + t270;
t68 = (-pkin(11) * qJD(4) - qJ(6) * qJD(5)) * t277 + (-qJD(6) * t178 + (-pkin(11) * qJD(5) - qJ(6) * qJD(4)) * t182) * t177 - t270;
t67 = -qJD(5) * t243 + t139 * t254 + t177 * t187 - t181 * t233;
t65 = -t178 * t252 + (pkin(5) * t178 - qJ(6) * t272) * qJD(4) + (-t163 + (-t151 + t283) * t177) * qJD(5) + t269;
t62 = pkin(5) * t105 + t82;
t40 = pkin(5) * t67 + t57;
t37 = (t103 * t175 + t123 * t173) * t183 + t298;
t21 = qJD(5) * t59 + t44 * t177 - t181 * t73;
t18 = pkin(5) * t58 + t24;
t17 = -qJ(6) * t67 - qJD(6) * t105 - t19;
t16 = t104 * pkin(5) + t66 * qJ(6) - t106 * qJD(6) + t20;
t5 = pkin(5) * t21 + t7;
t2 = -qJ(6) * t21 - qJD(6) * t58 - t3;
t1 = t43 * pkin(5) - t22 * qJ(6) - t59 * qJD(6) + t4;
t6 = [0, 0, 0, 0.2e1 * t180 * t242, 0.2e1 * (-t180 ^ 2 + t184 ^ 2) * t168 * qJD(2), 0.2e1 * t176 * t230, -0.2e1 * t176 * t231, 0, -0.2e1 * pkin(1) * t168 * t265 + 0.2e1 * t134 * t176, -0.2e1 * pkin(1) * t242 + 0.2e1 * t133 * t176, 0.2e1 * t99 * t74, -0.2e1 * t73 * t99 - 0.2e1 * t74 * t98, -0.2e1 * t137 * t74 + 0.2e1 * t221 * t99, 0.2e1 * t137 * t73 - 0.2e1 * t221 * t98, -0.2e1 * t137 * t221, -0.2e1 * t37 * t137 + 0.2e1 * t221 * t293 + 0.2e1 * t70 * t73 + 0.2e1 * t75 * t98, -0.2e1 * t137 * t36 - 0.2e1 * t221 * t247 + 0.2e1 * t70 * t74 + 0.2e1 * t75 * t99, 0.2e1 * t80 * t44, -0.2e1 * t43 * t80 - 0.2e1 * t44 * t79, 0.2e1 * t44 * t98 + 0.2e1 * t73 * t80, -0.2e1 * t43 * t98 - 0.2e1 * t73 * t79, 0.2e1 * t98 * t73, 0.2e1 * t10 * t98 + 0.2e1 * t26 * t73 + 0.2e1 * t35 * t79 + 0.2e1 * t43 * t52, -0.2e1 * t27 * t73 + 0.2e1 * t35 * t80 + 0.2e1 * t44 * t52 + 0.2e1 * t9 * t98, 0.2e1 * t59 * t22, -0.2e1 * t21 * t59 - 0.2e1 * t22 * t58, 0.2e1 * t22 * t79 + 0.2e1 * t43 * t59, -0.2e1 * t21 * t79 - 0.2e1 * t43 * t58, 0.2e1 * t79 * t43, 0.2e1 * t12 * t43 + 0.2e1 * t21 * t24 + 0.2e1 * t4 * t79 + 0.2e1 * t58 * t7, -0.2e1 * t13 * t43 + 0.2e1 * t22 * t24 + 0.2e1 * t3 * t79 + 0.2e1 * t59 * t7, -0.2e1 * t1 * t59 - 0.2e1 * t11 * t21 - 0.2e1 * t2 * t58 - 0.2e1 * t22 * t8, 0.2e1 * t1 * t8 + 0.2e1 * t11 * t2 + 0.2e1 * t18 * t5; 0, 0, 0, 0, 0, t230, -t231, 0, t134, t133 (t179 * t74 + t261 * t99) * t173 (-t179 * t73 + t183 * t74 + (-t179 * t99 - t183 * t98) * qJD(3)) * t173, -t137 * t232 + t175 * t74 + t179 * t223, t137 * t233 - t175 * t73 + t183 * t223, t175 * t221, t132 * t137 + t37 * t175 + ((pkin(2) * t278 - t159) * t231 - pkin(2) * t73 - t75 * t183 + t70 * t263) * t173, -t131 * t137 + t36 * t175 + (-pkin(2) * t74 + t179 * t75 - t231 * t268 + t261 * t70) * t173, t44 * t139 + t187 * t80, -t80 * t104 - t44 * t138 - t139 * t43 - t187 * t79, t98 * t239 + t139 * t73 + (-t44 * t183 + t195 * t98 + t263 * t80) * t173, -t104 * t98 - t138 * t73 + (t183 * t43 - t263 * t79) * t173 (-t183 * t73 + t263 * t98) * t173, t52 * t104 + t127 * t43 + t132 * t79 + t35 * t138 + t61 * t98 + t86 * t73 + (-t10 * t183 + t26 * t263) * t173, t52 * t239 + t127 * t44 + t132 * t80 + t35 * t139 + t60 * t98 - t87 * t73 + (-t9 * t183 + t195 * t52 - t263 * t27) * t173, t106 * t22 - t59 * t66, -t105 * t22 - t106 * t21 + t58 * t66 - t59 * t67, t104 * t59 + t106 * t43 + t138 * t22 - t66 * t79, -t104 * t58 - t105 * t43 - t138 * t21 - t67 * t79, t104 * t79 + t138 * t43, t104 * t12 + t105 * t7 + t138 * t4 + t20 * t79 + t21 * t82 + t24 * t67 + t43 * t54 + t57 * t58, -t104 * t13 + t106 * t7 + t138 * t3 + t19 * t79 + t22 * t82 - t24 * t66 - t43 * t55 + t57 * t59, -t1 * t106 - t105 * t2 - t11 * t67 - t16 * t59 - t17 * t58 - t21 * t42 - t22 * t41 + t66 * t8, t1 * t41 + t11 * t17 + t16 * t8 + t18 * t40 + t2 * t42 + t5 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t222, 0.2e1 * (-t179 ^ 2 + t183 ^ 2) * t167 * qJD(3), 0.2e1 * t175 * t232, -0.2e1 * t175 * t233, 0, -0.2e1 * pkin(2) * t167 * t263 - 0.2e1 * t132 * t175, -0.2e1 * pkin(2) * t241 + 0.2e1 * t131 * t175, 0.2e1 * t139 * t187, -0.2e1 * t139 * t104 - 0.2e1 * t138 * t187 (t139 * t263 - t175 * t236 - t195 * t281) * t292 (t104 * t183 - t138 * t263) * t292, -0.2e1 * t222, 0.2e1 * t127 * t104 + 0.2e1 * t132 * t138 + 0.2e1 * (-t183 * t61 + t263 * t86) * t173, 0.2e1 * t127 * t239 + 0.2e1 * t132 * t139 + 0.2e1 * (t127 * t195 - t60 * t183 - t263 * t87) * t173, -0.2e1 * t106 * t66, 0.2e1 * t105 * t66 - 0.2e1 * t106 * t67, 0.2e1 * t104 * t106 - 0.2e1 * t138 * t66, -0.2e1 * t104 * t105 - 0.2e1 * t138 * t67, 0.2e1 * t138 * t104, 0.2e1 * t104 * t54 + 0.2e1 * t105 * t57 + 0.2e1 * t138 * t20 + 0.2e1 * t67 * t82, -0.2e1 * t104 * t55 + 0.2e1 * t106 * t57 + 0.2e1 * t138 * t19 - 0.2e1 * t66 * t82, -0.2e1 * t105 * t17 - 0.2e1 * t106 * t16 + 0.2e1 * t41 * t66 - 0.2e1 * t42 * t67, 0.2e1 * t16 * t41 + 0.2e1 * t17 * t42 + 0.2e1 * t40 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t73, t221, t37, t36, t178 * t44 + t257 * t80, -t178 * t43 + t44 * t182 + (-t178 * t80 - t182 * t79) * qJD(4), t206, -t205, 0, -pkin(3) * t43 - pkin(11) * t206 - t35 * t182 + t259 * t52, -pkin(3) * t44 + pkin(11) * t205 + t35 * t178 + t257 * t52, t59 * t237 + (-t255 * t59 + t285) * t178, t212 * t257 + (-t286 - t181 * t21 + (t177 * t58 - t181 * t59) * qJD(5)) * t178 (t258 * t79 - t22) * t182 + (qJD(4) * t59 - t203) * t178 (-t260 * t79 + t21) * t182 + (-qJD(4) * t58 - t204) * t178, -t182 * t43 + t259 * t79, t124 * t43 + t85 * t79 + (-t4 + (pkin(11) * t58 + t177 * t24) * qJD(4)) * t182 + (pkin(11) * t21 + qJD(4) * t12 + t208) * t178, -t125 * t43 + t84 * t79 + (-t3 + (pkin(11) * t59 + t181 * t24) * qJD(4)) * t182 + (pkin(11) * t22 - qJD(4) * t13 - t207) * t178, -t108 * t21 - t97 * t22 - t68 * t58 - t65 * t59 + t214 * t257 + (-t1 * t181 - t177 * t2 + (-t11 * t181 + t177 * t8) * qJD(5)) * t178, t1 * t97 + t108 * t2 + t11 * t68 + t118 * t18 + t147 * t5 + t65 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, -t233, 0, -t132, t131, t178 * t219 + (-t138 * t178 + t139 * t182) * qJD(4), t172 * t232 - t178 * t104 + (-0.2e1 * t182 * t138 - t139 * t178) * qJD(4), t296 (t178 * t256 + t179 * t262) * t173, 0, -pkin(3) * t104 - pkin(11) * t296 + t127 * t259 - t132 * t182, t132 * t178 + t217 * t173 * t262 + (-pkin(11) * t178 * t281 + pkin(3) * t138 + t127 * t182) * qJD(4), t106 * t193 - t277 * t66, t210 * t257 + (t284 - t181 * t67 + (t105 * t177 - t106 * t181) * qJD(5)) * t178 (t138 * t258 + t66) * t182 + (qJD(4) * t106 - t197) * t178 (-t138 * t260 + t67) * t182 + (-qJD(4) * t105 - t198) * t178, -t104 * t182 + t138 * t259, t124 * t104 + t85 * t138 + (-t20 + (pkin(11) * t105 + t177 * t82) * qJD(4)) * t182 + (pkin(11) * t67 + qJD(4) * t54 + t202) * t178, -t125 * t104 + t84 * t138 + (-t19 + (pkin(11) * t106 + t181 * t82) * qJD(4)) * t182 + (-pkin(11) * t66 - qJD(4) * t55 - t201) * t178, -t68 * t105 - t65 * t106 - t108 * t67 + t97 * t66 + t213 * t257 + (-t16 * t181 - t17 * t177 + (t177 * t41 - t181 * t42) * qJD(5)) * t178, t108 * t17 + t118 * t62 + t147 * t40 + t16 * t97 + t41 * t65 + t42 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, -0.2e1 * t225, 0, 0, 0, t178 * t251, t182 * t251, -0.2e1 * t170 * t229 + 0.2e1 * t171 * t228, 0.2e1 * t170 * t226 + t220 * t297, 0.2e1 * t178 * t235 + 0.2e1 * t258 * t266, -0.2e1 * t177 * t225 + 0.2e1 * t178 * t234, -0.2e1 * t228, 0.2e1 * t124 * t259 - 0.2e1 * t85 * t182 + 0.2e1 * (t170 * t254 + t177 * t224) * pkin(11), -0.2e1 * t125 * t259 - 0.2e1 * t84 * t182 + 0.2e1 * (-t170 * t255 + t181 * t224) * pkin(11), 0.2e1 * (-t108 * t177 - t181 * t97) * t257 + 0.2e1 * (-t177 * t68 - t181 * t65 + (-t108 * t181 + t177 * t97) * qJD(5)) * t178, 0.2e1 * t108 * t68 + 0.2e1 * t118 * t147 + 0.2e1 * t65 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, t73, t10, t9, t254 * t59 + t286, qJD(5) * t212 - t177 * t21 + t285, t204, -t203, 0, -pkin(4) * t21 - pkin(12) * t204 + t207, -pkin(4) * t22 + pkin(12) * t203 + t208, qJD(5) * t214 - t1 * t177 - t135 * t58 - t136 * t59 - t152 * t22 + t153 * t21 + t2 * t181, t1 * t152 + t11 * t135 + t136 * t8 - t153 * t2 + t166 * t5 + t18 * t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, -t104, t233, t61, t60, t106 * t254 - t284, qJD(5) * t210 - t177 * t67 - t66 * t181, t198, -t197, 0, -pkin(4) * t67 - pkin(12) * t198 + t201, pkin(4) * t66 + pkin(12) * t197 + t202, qJD(5) * t213 - t135 * t105 - t136 * t106 + t152 * t66 + t153 * t67 - t16 * t177 + t17 * t181, t135 * t42 + t136 * t41 + t152 * t16 - t153 * t17 + t166 * t40 + t245 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, -t259, 0, -t246, pkin(11) * t259, -t178 * t226 + t220, t229 * t297 - t257 * t267, -t234 + t238, t192, 0 (pkin(12) * t272 + (-pkin(4) * t181 + t289) * t178) * qJD(5) + (t177 * t216 - t163) * qJD(4) (pkin(11) * t277 + t177 * t215) * qJD(5) + (t181 * t216 + t248) * qJD(4) (-t152 * t257 - t136 * t178 + t68 + (t153 * t178 - t97) * qJD(5)) * t181 + (t153 * t257 - t135 * t178 - t65 + (t152 * t178 - t108) * qJD(5)) * t177, t108 * t135 + t118 * t166 + t136 * t97 + t147 * t245 + t152 * t65 - t153 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t229, -0.2e1 * t226, 0, 0, 0, t177 * t250, t181 * t250, 0.2e1 * t135 * t181 - 0.2e1 * t136 * t177 + 0.2e1 * (-t152 * t181 + t153 * t177) * qJD(5), -0.2e1 * t135 * t153 + 0.2e1 * t136 * t152 + 0.2e1 * t166 * t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t43, t4, t3, -pkin(5) * t22, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t67, t104, t20, t19, pkin(5) * t66, t16 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t191, t259, t85, t84, -t193 * pkin(5), t65 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, -t255, 0, -pkin(12) * t254, pkin(12) * t255, -pkin(5) * t254, t136 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
