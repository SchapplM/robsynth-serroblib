% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRR6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:19
% EndTime: 2019-12-05 17:10:25
% DurationCPUTime: 2.36s
% Computational Cost: add. (3472->307), mult. (5904->419), div. (0->0), fcn. (4272->14), ass. (0->199)
t150 = qJD(2) + qJD(3);
t157 = sin(qJ(5));
t161 = cos(qJ(5));
t148 = qJDD(2) + qJDD(3);
t162 = cos(qJ(4));
t242 = t162 * t148;
t158 = sin(qJ(4));
t245 = t158 * t148;
t149 = qJD(4) + qJD(5);
t244 = t161 * t158;
t92 = t157 * t162 + t244;
t283 = t149 * t92;
t30 = t150 * t283 + t157 * t245 - t161 * t242;
t151 = t158 ^ 2;
t152 = t162 ^ 2;
t238 = t151 + t152;
t164 = cos(qJ(2));
t230 = t164 * qJD(1);
t125 = qJD(2) * pkin(2) + t230;
t159 = sin(qJ(3));
t163 = cos(qJ(3));
t235 = qJD(3) * t163;
t160 = sin(qJ(2));
t229 = qJD(1) * qJD(2);
t279 = t160 * qJDD(1) + t164 * t229;
t139 = t164 * qJDD(1);
t237 = qJD(1) * t160;
t282 = -qJDD(2) * pkin(2) + qJD(3) * t237 + t160 * t229 - t139;
t194 = -t125 * t235 + t282 * t159 - t279 * t163;
t33 = t148 * pkin(7) - t194;
t208 = t238 * t33;
t154 = qJ(2) + qJ(3);
t143 = cos(t154);
t270 = g(3) * t143;
t236 = qJD(3) * t159;
t193 = t125 * t236 + t279 * t159 + t282 * t163;
t269 = t148 * pkin(3);
t34 = t193 - t269;
t281 = t34 + t270;
t141 = sin(t154);
t156 = cos(pkin(9));
t254 = t141 * t156;
t155 = sin(pkin(9));
t255 = t141 * t155;
t280 = g(1) * t254 + g(2) * t255;
t231 = qJD(5) * t161;
t233 = qJD(4) * t162;
t278 = -t161 * t233 - t162 * t231;
t277 = -t270 + t280;
t243 = t161 * t162;
t246 = t157 * t158;
t90 = -t243 + t246;
t93 = t159 * t164 + t163 * t160;
t49 = t90 * t93;
t91 = t159 * t160 - t163 * t164;
t276 = t150 * t91;
t225 = pkin(2) * t235;
t128 = t159 * t237;
t82 = t163 * t230 - t128;
t275 = t225 - t82;
t81 = t93 * qJD(1);
t202 = pkin(2) * t236 - t81;
t199 = g(1) * t156 + g(2) * t155;
t136 = t159 * pkin(2) + pkin(7);
t267 = t163 * pkin(2);
t138 = -pkin(3) - t267;
t166 = qJD(4) ^ 2;
t274 = t136 * t166 + t138 * t148 + t150 * t202;
t218 = t150 * t233;
t73 = t159 * t125 + t163 * t237;
t67 = t150 * pkin(7) + t73;
t11 = -t67 * t233 + qJDD(4) * pkin(4) - t158 * t33 + (-t218 - t245) * pkin(8);
t234 = qJD(4) * t158;
t185 = t150 * t234 - t242;
t14 = -pkin(8) * t185 + t162 * t33 - t67 * t234;
t232 = qJD(5) * t157;
t215 = pkin(8) * t150 + t67;
t50 = t215 * t158;
t47 = qJD(4) * pkin(4) - t50;
t51 = t215 * t162;
t3 = (qJD(5) * t47 + t14) * t161 + t157 * t11 - t51 * t232;
t165 = -pkin(8) - pkin(7);
t273 = pkin(2) * t160;
t133 = g(3) * t141;
t268 = t150 * pkin(3);
t222 = t150 * t246;
t74 = -t150 * t243 + t222;
t76 = t92 * t150;
t266 = t76 * t74;
t265 = -pkin(8) - t136;
t84 = t265 * t158;
t144 = t162 * pkin(8);
t85 = t162 * t136 + t144;
t52 = -t157 * t85 + t161 * t84;
t210 = qJD(4) * t265;
t64 = t158 * t210 + t162 * t225;
t65 = -t158 * t225 + t162 * t210;
t264 = qJD(5) * t52 + t157 * t65 + t161 * t64 + t90 * t82;
t53 = t157 * t84 + t161 * t85;
t263 = -qJD(5) * t53 - t157 * t64 + t161 * t65 + t92 * t82;
t112 = t165 * t158;
t113 = t162 * pkin(7) + t144;
t62 = t161 * t112 - t157 * t113;
t72 = t163 * t125 - t128;
t219 = qJD(4) * t165;
t94 = t158 * t219;
t95 = t162 * t219;
t262 = qJD(5) * t62 + t157 * t95 + t161 * t94 + t90 * t72;
t63 = t157 * t112 + t161 * t113;
t261 = -qJD(5) * t63 - t157 * t94 + t161 * t95 + t92 * t72;
t223 = pkin(4) * t234;
t260 = t223 + t202;
t259 = t157 * t51;
t258 = t161 * t51;
t257 = t276 * t150;
t256 = t73 * t150;
t253 = t143 * t155;
t252 = t143 * t156;
t251 = t150 * t158;
t153 = qJ(4) + qJ(5);
t140 = sin(t153);
t250 = t155 * t140;
t142 = cos(t153);
t249 = t155 * t142;
t248 = t156 * t140;
t247 = t156 * t142;
t241 = qJDD(1) - g(3);
t240 = t143 * pkin(3) + t141 * pkin(7);
t239 = t151 - t152;
t66 = -t72 - t268;
t227 = t280 * t162 + t66 * t234;
t226 = t281 * t158 + t66 * t233;
t146 = t150 ^ 2;
t221 = t158 * t146 * t162;
t220 = g(1) * t252 + g(2) * t253 + t133;
t137 = t162 * pkin(4) + pkin(3);
t209 = t72 * t238;
t206 = t143 * t137 - t141 * t165;
t205 = t238 * t148;
t204 = -t148 * t244 + t278 * t150 - t157 * t242;
t203 = t158 * t218;
t201 = -t73 + t223;
t200 = -pkin(3) * t141 - t273;
t198 = g(1) * t155 - g(2) * t156;
t195 = -t220 + t208;
t192 = t149 * t246;
t58 = t150 * t93;
t191 = -t91 * t148 - t58 * t150;
t22 = t157 * t47 + t258;
t189 = t137 * t141 + t143 * t165;
t188 = t199 * t141;
t187 = t198 * t162;
t21 = t161 * t47 - t259;
t4 = -qJD(5) * t22 + t161 * t11 - t157 * t14;
t55 = t192 + t278;
t186 = t21 * t55 - t22 * t283 - t3 * t90 - t4 * t92 - t220;
t25 = pkin(4) * t185 + t34;
t54 = -t137 * t150 - t72;
t184 = t277 * t142 + t25 * t90 + t283 * t54;
t183 = pkin(7) * t166 - t256 - t269;
t182 = t166 * t93 - t191;
t180 = t194 + t220;
t179 = -pkin(7) * qJDD(4) + (t72 - t268) * qJD(4);
t178 = 0.2e1 * t276 * qJD(4) - qJDD(4) * t93;
t177 = t199 * t143 + t133;
t176 = -g(3) * t164 + t199 * t160;
t175 = -t193 + t277;
t174 = t275 * t238;
t173 = t25 * t92 - t54 * t55 + (-t188 + t270) * t140;
t172 = -qJDD(4) * t136 + (t138 * t150 - t275) * qJD(4);
t171 = -t66 * t150 + t177 - t33;
t169 = -g(1) * (-t143 * t247 - t250) - g(2) * (-t143 * t249 + t248) + t54 * t74 + t142 * t133 - t3;
t168 = -g(1) * (-t143 * t248 + t249) - g(2) * (-t143 * t250 - t247) - t54 * t76 + t140 * t133 + t4;
t167 = qJD(2) ^ 2;
t147 = qJDD(4) + qJDD(5);
t145 = t164 * pkin(2);
t115 = pkin(7) * t252;
t114 = pkin(7) * t253;
t107 = -t137 - t267;
t106 = qJDD(4) * t162 - t166 * t158;
t105 = qJDD(4) * t158 + t166 * t162;
t78 = t152 * t148 - 0.2e1 * t203;
t77 = t151 * t148 + 0.2e1 * t203;
t61 = -0.2e1 * t239 * t150 * qJD(4) + 0.2e1 * t158 * t242;
t48 = t92 * t93;
t43 = -t90 * t147 - t149 * t283;
t42 = t92 * t147 - t55 * t149;
t35 = -t74 ^ 2 + t76 ^ 2;
t29 = t150 * t192 + t204;
t24 = -t161 * t50 - t259;
t23 = t157 * t50 - t258;
t20 = t76 * t149 - t30;
t19 = -t204 + (-t222 + t74) * t149;
t9 = t283 * t74 + t30 * t90;
t8 = -t29 * t92 - t76 * t55;
t7 = t149 * t49 + t276 * t92;
t6 = t276 * t90 - t283 * t93;
t5 = -t283 * t76 + t29 * t90 - t92 * t30 + t55 * t74;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t241, 0, 0, 0, 0, 0, 0, t164 * qJDD(2) - t167 * t160, -qJDD(2) * t160 - t167 * t164, 0, -g(3) + (t160 ^ 2 + t164 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t191, -t93 * t148 + t257, 0, t193 * t91 - t194 * t93 - t276 * t73 - t72 * t58 - g(3), 0, 0, 0, 0, 0, 0, t158 * t178 - t162 * t182, t158 * t182 + t162 * t178, t93 * t205 - t238 * t257, t34 * t91 + t66 * t58 - g(3) + t238 * (-t276 * t67 + t33 * t93), 0, 0, 0, 0, 0, 0, -t48 * t147 + t7 * t149 + t91 * t30 + t58 * t74, t49 * t147 - t6 * t149 - t91 * t29 + t58 * t76, -t48 * t29 + t49 * t30 - t6 * t74 - t7 * t76, t21 * t7 + t22 * t6 + t25 * t91 - t3 * t49 - t4 * t48 + t54 * t58 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t139 + t176, -t241 * t160 + t199 * t164, 0, 0, 0, 0, 0, 0, 0, t148, t81 * t150 + (t148 * t163 - t150 * t236) * pkin(2) + t175, t82 * t150 + (-t148 * t159 - t150 * t235) * pkin(2) + t180, 0, t72 * t81 - t73 * t82 + (-t159 * t194 - t163 * t193 + (-t159 * t72 + t163 * t73) * qJD(3) + t176) * pkin(2), t77, t61, t105, t78, t106, 0, t172 * t158 + (-t281 - t274) * t162 + t227, t172 * t162 + (-t188 + t274) * t158 + t226, t136 * t205 + t150 * t174 + t195, t34 * t138 - g(1) * (t156 * t200 + t115) - g(2) * (t155 * t200 + t114) - g(3) * (t145 + t240) + t202 * t66 + t136 * t208 + t174 * t67, t8, t5, t42, t9, t43, 0, t107 * t30 + t52 * t147 + t263 * t149 + t260 * t74 + t184, -t107 * t29 - t53 * t147 - t264 * t149 + t260 * t76 + t173, -t263 * t76 - t264 * t74 + t52 * t29 - t53 * t30 + t186, t3 * t53 + t4 * t52 + t25 * t107 - g(3) * (t145 + t206) + t260 * t54 + t264 * t22 + t263 * t21 + t199 * (t189 + t273); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t175 + t256, t72 * t150 + t180, 0, 0, t77, t61, t105, t78, t106, 0, t179 * t158 + (-t183 - t281) * t162 + t227, t179 * t162 + (-t188 + t183) * t158 + t226, pkin(7) * t205 - t150 * t209 + t195, -t34 * pkin(3) - t66 * t73 - g(1) * (-pkin(3) * t254 + t115) - g(2) * (-pkin(3) * t255 + t114) - g(3) * t240 - t67 * t209 + pkin(7) * t208, t8, t5, t42, t9, t43, 0, -t137 * t30 + t62 * t147 + t261 * t149 + t201 * t74 + t184, t137 * t29 - t63 * t147 - t262 * t149 + t201 * t76 + t173, -t261 * t76 - t262 * t74 + t62 * t29 - t63 * t30 + t186, -g(3) * t206 - t25 * t137 + t199 * t189 + t201 * t54 + t261 * t21 + t262 * t22 + t3 * t63 + t4 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221, t239 * t146, t245, t221, t242, qJDD(4), t158 * t171 - t187, t158 * t198 + t162 * t171, 0, 0, t266, t35, t19, -t266, t20, t147, -t23 * t149 + (t147 * t161 - t149 * t232 - t74 * t251) * pkin(4) + t168, t24 * t149 + (-t147 * t157 - t149 * t231 - t76 * t251) * pkin(4) + t169, (t22 + t23) * t76 + (-t21 + t24) * t74 + (-t157 * t30 + t161 * t29 + (t157 * t76 - t161 * t74) * qJD(5)) * pkin(4), -t21 * t23 - t22 * t24 + (t3 * t157 + t4 * t161 - t187 + (-t54 * t150 + t177) * t158 + (-t21 * t157 + t22 * t161) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, t35, t19, -t266, t20, t147, t22 * t149 + t168, t21 * t149 + t169, 0, 0;];
tau_reg = t1;
