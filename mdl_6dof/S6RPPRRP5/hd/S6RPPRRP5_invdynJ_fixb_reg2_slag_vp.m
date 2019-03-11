% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:53
% EndTime: 2019-03-09 02:08:59
% DurationCPUTime: 3.27s
% Computational Cost: add. (3432->436), mult. (6104->503), div. (0->0), fcn. (3454->6), ass. (0->231)
t134 = sin(qJ(5));
t137 = cos(qJ(5));
t223 = t137 * qJD(4);
t138 = cos(qJ(4));
t233 = qJD(1) * t138;
t79 = t134 * t233 - t223;
t135 = sin(qJ(4));
t234 = qJD(1) * t135;
t99 = qJD(5) + t234;
t166 = t79 * t99;
t224 = qJD(5) * t138;
t200 = t134 * t224;
t216 = t138 * qJDD(1);
t39 = (t135 * t223 + t200) * qJD(1) - qJD(5) * t223 - t134 * qJDD(4) - t137 * t216;
t289 = -t39 - t166;
t27 = -t39 + t166;
t133 = pkin(1) + qJ(3);
t179 = pkin(4) * t135 - pkin(8) * t138;
t87 = t179 + t133;
t57 = t87 * qJD(1) - qJD(2);
t108 = qJ(2) * qJD(1) + qJD(3);
t95 = -pkin(7) * qJD(1) + t108;
t86 = t135 * t95;
t63 = qJD(4) * pkin(8) + t86;
t32 = -t134 * t63 + t137 * t57;
t230 = qJD(4) * t134;
t81 = t137 * t233 + t230;
t29 = -qJ(6) * t81 + t32;
t20 = pkin(5) * t99 + t29;
t33 = t134 * t57 + t137 * t63;
t30 = -qJ(6) * t79 + t33;
t174 = t134 * t30 + t137 * t20;
t136 = sin(qJ(1));
t139 = cos(qJ(1));
t237 = g(1) * t139 + g(2) * t136;
t288 = -t174 * qJD(1) - t237;
t268 = t81 * t99;
t222 = qJD(1) * qJD(4);
t198 = t135 * t222;
t250 = qJD(5) * t81;
t40 = -t137 * qJDD(4) + (-t198 + t216) * t134 + t250;
t287 = t40 - t268;
t286 = t40 + t268;
t132 = -pkin(7) + qJ(2);
t228 = qJD(4) * t138;
t285 = qJD(2) * t135 + t132 * t228;
t128 = t135 ^ 2;
t129 = t138 ^ 2;
t236 = t128 + t129;
t284 = t236 * t95;
t283 = qJD(1) * t133;
t247 = t134 * t138;
t240 = t139 * t137;
t244 = t136 * t134;
t65 = t135 * t244 - t240;
t243 = t136 * t137;
t245 = t135 * t139;
t67 = -t134 * t245 - t243;
t282 = -g(1) * t67 + g(2) * t65 + g(3) * t247;
t96 = -qJD(2) + t283;
t281 = (qJD(2) + t96 + t283) * qJD(4) + qJDD(4) * t132;
t280 = t81 ^ 2;
t125 = qJD(1) * qJD(2);
t117 = 0.2e1 * t125;
t197 = t138 * t222;
t217 = t135 * qJDD(1);
t76 = qJDD(5) + t197 + t217;
t279 = pkin(5) * t76;
t278 = pkin(5) * t79;
t275 = pkin(5) * t134;
t274 = pkin(7) * t139;
t122 = g(1) * t136;
t121 = g(2) * t139;
t273 = g(3) * t135;
t272 = g(3) * t138;
t271 = t32 * t99;
t270 = t33 * t99;
t269 = t81 * t79;
t267 = qJ(6) + pkin(8);
t266 = -t29 + t20;
t192 = qJD(5) * t267;
t203 = t134 * t234;
t242 = t137 * t138;
t180 = pkin(4) * t138 + pkin(8) * t135;
t85 = t180 * qJD(1);
t48 = t134 * t85 + t95 * t242;
t265 = -qJ(6) * t203 + qJD(6) * t137 - t134 * t192 - t48;
t246 = t135 * t137;
t47 = t137 * t85 - t95 * t247;
t264 = -qJD(6) * t134 - t137 * t192 - (pkin(5) * t138 + qJ(6) * t246) * qJD(1) - t47;
t241 = t138 * t139;
t263 = g(1) * t134 * t241 + g(2) * t138 * t244;
t53 = t132 * t246 + t134 * t87;
t262 = qJ(6) * t39;
t261 = qJ(6) * t40;
t260 = t134 * t39;
t259 = t137 * t40;
t258 = t137 * t81;
t257 = t138 * t79;
t256 = t138 * t81;
t126 = qJ(2) * qJDD(1);
t194 = qJDD(3) + t125 + t126;
t83 = -pkin(7) * qJDD(1) + t194;
t255 = t138 * t83;
t254 = qJD(1) * t96;
t253 = qJD(1) * t99;
t252 = qJD(4) * t79;
t251 = qJD(4) * t81;
t141 = qJD(1) ^ 2;
t249 = t128 * t141;
t248 = t134 * t135;
t239 = t139 * pkin(1) + t136 * qJ(2);
t238 = t122 - t121;
t140 = qJD(4) ^ 2;
t235 = -t140 - t141;
t231 = qJD(2) * t138;
t229 = qJD(4) * t135;
t227 = qJD(5) * t132;
t226 = qJD(5) * t134;
t225 = qJD(5) * t137;
t221 = qJD(3) * qJD(1);
t220 = qJDD(1) * t133;
t218 = qJDD(4) * t135;
t130 = qJDD(1) * pkin(1);
t215 = t130 - qJDD(2);
t214 = pkin(8) * qJD(5) * t99;
t213 = t99 * t248;
t212 = t99 * t246;
t210 = t99 * t230;
t209 = t99 * t223;
t193 = -qJD(6) - t278;
t64 = -qJD(4) * pkin(4) - t138 * t95;
t45 = -t193 + t64;
t208 = t45 * t225;
t207 = t138 * t141 * t135;
t206 = t99 * t233;
t205 = t139 * qJ(3) + t239;
t204 = pkin(7) + t275;
t201 = t135 * t227;
t199 = t137 * t224;
t196 = qJDD(2) - t238;
t124 = qJDD(1) * qJ(3);
t195 = -t124 - t215;
t191 = t236 * t83;
t77 = t180 * qJD(4) + qJD(3);
t38 = t77 * qJD(1) + t179 * qJDD(1) - t195;
t44 = qJDD(4) * pkin(8) + t135 * t83 + t95 * t228;
t7 = t134 * t38 + t137 * t44 + t57 * t225 - t63 * t226;
t190 = t134 * t77 + t285 * t137 + t87 * t225;
t189 = -qJDD(4) * pkin(4) + t95 * t229;
t188 = qJD(5) * t79 - t39;
t187 = -t40 + t250;
t186 = 0.2e1 * t126 + t117 - t237;
t185 = g(2) * t205;
t184 = t135 * t197;
t183 = -t130 + t196;
t182 = -g(1) * t65 - g(2) * t67;
t66 = -t134 * t139 - t135 * t243;
t68 = t135 * t240 - t244;
t181 = -g(1) * t66 - g(2) * t68;
t177 = -t237 - t254;
t116 = t139 * qJ(2);
t175 = -t133 * t136 + t116;
t173 = t134 * t20 - t137 * t30;
t172 = t134 * t33 + t137 * t32;
t171 = t134 * t32 - t137 * t33;
t84 = -t195 + t221;
t170 = t96 * qJD(3) + t84 * t133;
t107 = pkin(5) * t137 + pkin(4);
t168 = t107 * t135 - t138 * t267;
t165 = -t124 + t183;
t164 = t237 - t83;
t162 = t134 * t76 + t99 * t225;
t161 = -t137 * t76 + t99 * t226;
t160 = t237 * t138;
t156 = pkin(5) * t40 + qJDD(6) + t189;
t43 = t189 - t255;
t154 = -t160 - t43;
t153 = -pkin(8) * t76 + t99 * t64;
t152 = t164 + t254;
t151 = g(1) * t68 - g(2) * t66 + g(3) * t242 - t7;
t150 = -t135 * t237 - t272;
t149 = -qJD(6) * t138 + (qJ(6) * qJD(4) - t227) * t135;
t8 = -qJD(5) * t33 - t134 * t44 + t137 * t38;
t148 = -t132 * t140 + t122 + t220 + t221 + t84;
t2 = -qJD(6) * t81 + t262 + t279 + t8;
t3 = -qJD(6) * t79 - t261 + t7;
t147 = -t174 * qJD(5) - t134 * t2 + t137 * t3;
t146 = -t172 * qJD(5) - t8 * t134 + t7 * t137;
t145 = t173 * qJD(5) - t134 * t3 - t137 * t2;
t144 = t171 * qJD(5) - t134 * t7 - t137 * t8;
t142 = t8 + t282;
t113 = t129 * t141;
t111 = qJDD(4) * t138;
t103 = g(2) * t241;
t100 = g(3) * t246;
t90 = t267 * t137;
t89 = t267 * t134;
t75 = t79 ^ 2;
t74 = (-t132 + t275) * t138;
t72 = t137 * t87;
t61 = t137 * t77;
t55 = -pkin(5) * t203 + t86;
t52 = -t132 * t248 + t72;
t50 = t135 * t76 + t99 * t228;
t49 = t132 * t229 - t231 + (-t134 * t229 + t199) * pkin(5);
t46 = -qJ(6) * t247 + t53;
t42 = -qJ(6) * t242 + t72 + (-t132 * t134 + pkin(5)) * t135;
t36 = -t75 + t280;
t26 = (t212 + t256) * qJD(1) + t162;
t25 = (t212 - t256) * qJD(1) + t162;
t24 = (-t213 + t257) * qJD(1) - t161;
t23 = (t213 + t257) * qJD(1) + t161;
t22 = -t137 * t201 + t61 + (-qJD(5) * t87 - t285) * t134;
t21 = -t134 * t201 + t190;
t19 = t156 - t255;
t18 = t134 * t166 - t259;
t17 = t99 * t258 - t260;
t16 = t79 * t199 + (t138 * t40 - t79 * t229) * t134;
t15 = -t81 * t200 + (-t138 * t39 - t81 * t229) * t137;
t14 = -qJ(6) * t199 + t149 * t134 + t190;
t13 = pkin(5) * t228 + t61 + t149 * t137 + ((qJ(6) * t138 - t87) * qJD(5) - t285) * t134;
t12 = (-t40 + t210) * t135 + (-t162 - t252) * t138;
t11 = (-t39 - t209) * t135 + (-t161 + t251) * t138;
t10 = -t137 * t253 + (-t40 - t210) * t138 + (-t162 + t252) * t135;
t9 = t134 * t253 + (t39 - t209) * t138 + (t161 + t251) * t135;
t6 = t287 * t134 + t27 * t137;
t5 = -t286 * t134 + t289 * t137;
t4 = (t134 * t81 + t137 * t79) * t229 + (t260 - t259 + (t134 * t79 - t258) * qJD(5)) * t138;
t1 = (qJD(1) * t81 + t187 * t135 - t79 * t228) * t137 + (qJD(1) * t79 + t188 * t135 + t81 * t228) * t134;
t28 = [0, 0, 0, 0, 0, qJDD(1), t238, t237, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t130 + t196, t186, t215 * pkin(1) - g(1) * (-pkin(1) * t136 + t116) - g(2) * t239 + (t117 + t126) * qJ(2), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t186, -t165 + t220 + 0.2e1 * t221, -g(1) * t175 + t194 * qJ(2) + t108 * qJD(2) + t170 - t185, qJDD(1) * t129 - 0.2e1 * t184, -0.2e1 * t135 * t216 + 0.2e1 * (t128 - t129) * t222, -t135 * t140 + t111, qJDD(1) * t128 + 0.2e1 * t184, -t138 * t140 - t218, 0, t281 * t138 + (t148 - t121) * t135, -t135 * t281 + t148 * t138 - t103, t237 + t236 * (-qJDD(1) * t132 - t125 - t83) -g(1) * (t175 - t274) - g(2) * (-pkin(7) * t136 + t205) + t132 * t191 + qJD(2) * t284 + t170, t15, t4, t11, t16, t12, t50, t22 * t99 + t52 * t76 + (t8 + (t132 * t79 - t134 * t64) * qJD(4)) * t135 + (-qJD(2) * t79 + qJD(4) * t32 - t132 * t40 + t134 * t43 + t225 * t64) * t138 + t181, -t21 * t99 - t53 * t76 + (-t7 + (t132 * t81 - t137 * t64) * qJD(4)) * t135 + (-qJD(2) * t81 - qJD(4) * t33 + t132 * t39 + t137 * t43 - t226 * t64) * t138 + t182, -t21 * t79 - t22 * t81 + t39 * t52 - t40 * t53 + t103 + t172 * t229 + (t144 - t122) * t138, t7 * t53 + t33 * t21 + t8 * t52 + t32 * t22 - t64 * t231 - g(1) * (t116 - t274) - g(2) * (pkin(4) * t245 - pkin(8) * t241 + t205) + (-t138 * t43 + t229 * t64) * t132 + (g(2) * pkin(7) + g(1) * t87) * t136, t15, t4, t11, t16, t12, t50, t13 * t99 + t40 * t74 + t42 * t76 + t49 * t79 + (-t230 * t45 + t2) * t135 + (qJD(4) * t20 + t134 * t19 + t208) * t138 + t181, -t14 * t99 - t39 * t74 - t46 * t76 + t49 * t81 + (-t223 * t45 - t3) * t135 + (-qJD(4) * t30 + t137 * t19 - t226 * t45) * t138 + t182, -t13 * t81 - t14 * t79 + t39 * t42 - t40 * t46 + t103 + t174 * t229 + (t145 - t122) * t138, t3 * t46 + t30 * t14 + t2 * t42 + t20 * t13 + t19 * t74 + t45 * t49 - g(1) * t116 - t185 + (g(1) * t204 - g(2) * t168) * t139 + (-g(1) * (-t168 - t133) + g(2) * t204) * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t141, -qJ(2) * t141 + t183, 0, 0, 0, 0, 0, 0, 0, -t141, -qJDD(1) (-qJD(3) - t108) * qJD(1) + t165, 0, 0, 0, 0, 0, 0, -0.2e1 * t197 - t217, 0.2e1 * t198 - t216, t113 + t249 (-qJD(3) - t284) * qJD(1) + t165, 0, 0, 0, 0, 0, 0, t23, t26, t6 (t135 * t171 + t138 * t64) * qJD(1) + t144 - t238, 0, 0, 0, 0, 0, 0, t23, t26, t6 (t135 * t173 + t138 * t45) * qJD(1) + t145 - t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t141, t177 + t194, 0, 0, 0, 0, 0, 0, t235 * t135 + t111, t235 * t138 - t218, -t236 * qJDD(1), t191 + t177, 0, 0, 0, 0, 0, 0, t10, t9, t1, -t172 * qJD(1) + (-qJD(4) * t171 - t43) * t138 + (qJD(4) * t64 + t146) * t135 - t237, 0, 0, 0, 0, 0, 0, t10, t9, t1 (-qJD(4) * t173 - t19) * t138 + (qJD(4) * t45 + t147) * t135 + t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, t113 - t249, t216, -t207, -t217, qJDD(4), -t138 * t152 + t273, t135 * t152 + t272, 0, 0, t17, t5, t25, t18, t24, -t206, -t32 * t233 - t79 * t86 - pkin(4) * t40 - t47 * t99 + t100 + t153 * t134 + (t154 - t214) * t137, t33 * t233 - t81 * t86 + pkin(4) * t39 + t48 * t99 + t153 * t137 + (t214 + t43 - t273) * t134 + t263, t47 * t81 + t48 * t79 + (pkin(8) * t187 - t271 + t7) * t137 + (pkin(8) * t188 - t270 - t8) * t134 + t150, -t64 * t86 - t32 * t47 - t33 * t48 + (t154 + t273) * pkin(4) + (t146 + t150) * pkin(8), t17, t5, t25, t18, t24, -t206, -t20 * t233 - t107 * t40 - t55 * t79 - t76 * t89 + t100 + t264 * t99 + (-t19 - t160) * t137 + (t45 * t234 + (t45 + t278) * qJD(5)) * t134, t208 + t107 * t39 - t55 * t81 - t76 * t90 - t265 * t99 + (t138 * t30 + t246 * t45) * qJD(1) + (pkin(5) * t250 + t19 - t273) * t134 + t263, t288 * t135 - t264 * t81 - t265 * t79 - t39 * t89 - t40 * t90 + t147 - t272, t3 * t90 - t2 * t89 - t19 * t107 + g(3) * t168 + (pkin(5) * t226 - t55) * t45 + t265 * t30 + t264 * t20 - t237 * (t107 * t138 + t135 * t267); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t269, t36, t27, -t269, -t287, t76, -t64 * t81 + t142 + t270, t64 * t79 + t151 + t271, 0, 0, t269, t36, t27, -t269, -t287, t76, 0.2e1 * t279 + t262 + t30 * t99 + (t193 - t45) * t81 + t142, -pkin(5) * t280 + t261 + t29 * t99 + (qJD(6) + t45) * t79 + t151, pkin(5) * t39 - t266 * t79, t266 * t30 + (-t45 * t81 + t2 + t282) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t286, t289, -t75 - t280, t138 * t164 + t20 * t81 + t30 * t79 + t156 - t273;];
tau_reg  = t28;
