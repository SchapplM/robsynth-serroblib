% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:02
% EndTime: 2019-12-05 18:46:11
% DurationCPUTime: 3.86s
% Computational Cost: add. (16146->360), mult. (36103->465), div. (0->0), fcn. (25891->8), ass. (0->233)
t201 = sin(qJ(3));
t205 = cos(qJ(3));
t206 = cos(qJ(2));
t202 = sin(qJ(2));
t247 = qJD(1) * t202;
t172 = qJD(1) * t205 * t206 - t201 * t247;
t173 = (t201 * t206 + t202 * t205) * qJD(1);
t200 = sin(qJ(4));
t204 = cos(qJ(4));
t153 = -t172 * t204 + t173 * t200;
t155 = t172 * t200 + t173 * t204;
t121 = t155 * t153;
t196 = qJDD(2) + qJDD(3);
t188 = qJDD(4) + t196;
t289 = t121 - t188;
t293 = t289 * pkin(4);
t197 = qJD(2) + qJD(3);
t189 = qJD(4) + t197;
t141 = pkin(4) * t189 - qJ(5) * t155;
t190 = t202 * qJDD(1);
t242 = qJD(1) * qJD(2);
t232 = t206 * t242;
t179 = t190 + t232;
t191 = t206 * qJDD(1);
t233 = t202 * t242;
t180 = t191 - t233;
t228 = t201 * t179 - t180 * t205;
t137 = -qJD(3) * t173 - t228;
t199 = t206 ^ 2;
t208 = qJD(1) ^ 2;
t203 = sin(qJ(1));
t277 = cos(qJ(1));
t231 = t203 * g(1) - g(2) * t277;
t219 = qJDD(1) * pkin(1) + t231;
t220 = qJD(2) * pkin(2) - pkin(7) * t247;
t140 = t180 * pkin(2) - t220 * t247 + t219 + (pkin(7) * t199 + pkin(6)) * t208;
t170 = t172 ^ 2;
t224 = pkin(3) * t197 - pkin(8) * t173;
t90 = t137 * pkin(3) + t170 * pkin(8) - t173 * t224 + t140;
t209 = -t141 * t155 - qJDD(5) + t90;
t222 = t205 * t179 + t201 * t180;
t138 = qJD(3) * t172 + t222;
t230 = -t137 * t204 + t200 * t138;
t92 = -qJD(4) * t155 - t230;
t292 = -t92 * pkin(4) - t209;
t258 = t200 * t289;
t161 = t172 * t173;
t287 = t161 + t196;
t291 = t201 * t287;
t251 = t204 * t289;
t290 = t205 * t287;
t168 = t197 * t172;
t127 = t138 - t168;
t145 = t189 * t153;
t253 = t202 * t208;
t221 = g(1) * t277 + g(2) * t203;
t264 = qJDD(1) * pkin(6);
t175 = -pkin(1) * t208 - t221 + t264;
t255 = t202 * t175;
t131 = qJDD(2) * pkin(2) - t179 * pkin(7) - t255 + (pkin(2) * t253 + pkin(7) * t242 - g(3)) * t206;
t165 = -t202 * g(3) + t206 * t175;
t193 = t199 * t208;
t132 = -pkin(2) * t193 + t180 * pkin(7) - qJD(2) * t220 + t165;
t107 = -t131 * t205 + t201 * t132;
t70 = pkin(3) * t287 - pkin(8) * t127 - t107;
t108 = t201 * t131 + t205 * t132;
t72 = -t170 * pkin(3) + t137 * pkin(8) - t197 * t224 + t108;
t45 = t200 * t72 - t204 * t70;
t286 = qJ(5) * t145 + 0.2e1 * qJD(5) * t155 + t293 + t45;
t151 = t153 ^ 2;
t46 = t200 * t70 + t204 * t72;
t35 = -pkin(4) * t151 + qJ(5) * t92 - 0.2e1 * qJD(5) * t153 - t141 * t189 + t46;
t152 = t155 ^ 2;
t171 = t173 ^ 2;
t187 = t189 ^ 2;
t195 = t197 ^ 2;
t215 = (-qJD(4) + t189) * t155 - t230;
t223 = t200 * t137 + t204 * t138;
t93 = -qJD(4) * t153 + t223;
t84 = t145 + t93;
t52 = t200 * t215 - t204 * t84;
t54 = t200 * t84 + t204 * t215;
t31 = t201 * t54 + t205 * t52;
t284 = pkin(7) * t31;
t113 = -t187 - t151;
t88 = t113 * t200 - t251;
t89 = t113 * t204 + t258;
t57 = t201 * t89 + t205 * t88;
t283 = pkin(7) * t57;
t136 = -t152 - t187;
t115 = t121 + t188;
t259 = t200 * t115;
t98 = t136 * t204 - t259;
t252 = t204 * t115;
t99 = -t136 * t200 - t252;
t63 = t201 * t99 + t205 * t98;
t282 = pkin(7) * t63;
t281 = pkin(8) * t52;
t280 = pkin(8) * t88;
t279 = pkin(8) * t98;
t34 = -qJ(5) * t93 - t286;
t267 = t204 * t34;
t12 = t200 * t35 + t267;
t33 = pkin(4) * t34;
t276 = pkin(3) * t12 + t33;
t58 = -t201 * t88 + t205 * t89;
t243 = qJD(4) + t189;
t80 = t155 * t243 + t230;
t275 = pkin(6) * (-t202 * t57 + t206 * t58) - pkin(1) * t80;
t64 = -t201 * t98 + t205 * t99;
t85 = -t153 * t243 + t223;
t274 = pkin(6) * (-t202 * t63 + t206 * t64) - pkin(1) * t85;
t50 = pkin(3) * t52;
t77 = pkin(4) * t84;
t273 = -t77 + t50;
t106 = -t151 - t152;
t32 = -t201 * t52 + t205 * t54;
t272 = pkin(6) * (-t202 * t31 + t206 * t32) - pkin(1) * t106;
t271 = t200 * t34;
t270 = t200 * t90;
t25 = t200 * t46 - t204 * t45;
t269 = t201 * t25;
t66 = -t107 * t205 + t108 * t201;
t268 = t202 * t66;
t266 = t204 * t90;
t265 = t205 * t25;
t263 = t189 * t200;
t262 = t189 * t204;
t261 = t197 * t201;
t260 = t197 * t205;
t257 = t201 * t140;
t158 = -t161 + t196;
t256 = t201 * t158;
t185 = t206 * t253;
t254 = t202 * (qJDD(2) + t185);
t250 = t205 * t140;
t249 = t205 * t158;
t248 = t206 * (qJDD(2) - t185);
t244 = qJD(3) + t197;
t97 = pkin(3) * t98;
t241 = t97 - t46;
t240 = -pkin(2) * t80 + pkin(7) * t58;
t239 = -pkin(2) * t85 + pkin(7) * t64;
t238 = -pkin(3) * t80 + pkin(8) * t89;
t237 = -pkin(3) * t85 + pkin(8) * t99;
t235 = -pkin(2) * t106 + pkin(7) * t32;
t234 = -pkin(3) * t106 + pkin(8) * t54;
t26 = t200 * t45 + t204 * t46;
t67 = t107 * t201 + t108 * t205;
t164 = t206 * g(3) + t255;
t229 = t202 * t164 + t165 * t206;
t87 = pkin(3) * t88;
t227 = -t45 + t87;
t217 = pkin(4) * t136 - t35;
t216 = t97 + t217;
t214 = (-qJD(3) + t197) * t173 - t228;
t212 = t34 - t293;
t211 = t212 + t87;
t207 = qJD(2) ^ 2;
t198 = t202 ^ 2;
t192 = t198 * t208;
t181 = t191 - 0.2e1 * t233;
t178 = t190 + 0.2e1 * t232;
t174 = t208 * pkin(6) + t219;
t167 = -t171 + t195;
t166 = t170 - t195;
t163 = -t171 - t195;
t160 = t171 - t170;
t156 = -t195 - t170;
t143 = -t152 + t187;
t142 = t151 - t187;
t139 = -t170 - t171;
t129 = -t163 * t201 - t249;
t128 = t163 * t205 - t256;
t126 = t138 + t168;
t125 = t172 * t244 + t222;
t122 = t173 * t244 + t228;
t119 = t152 - t151;
t118 = t156 * t205 - t291;
t117 = t156 * t201 + t290;
t110 = (-t153 * t204 + t155 * t200) * t189;
t109 = (-t153 * t200 - t155 * t204) * t189;
t103 = t142 * t204 - t259;
t102 = -t143 * t200 - t251;
t101 = t142 * t200 + t252;
t100 = t143 * t204 - t258;
t95 = t127 * t201 + t205 * t214;
t94 = -t127 * t205 + t201 * t214;
t83 = -t145 + t93;
t76 = -t155 * t263 + t204 * t93;
t75 = t155 * t262 + t200 * t93;
t74 = t153 * t262 - t200 * t92;
t73 = t153 * t263 + t204 * t92;
t65 = -pkin(4) * t85 - qJ(5) * t115;
t62 = pkin(2) * t63;
t60 = -t266 - t279;
t59 = -t270 - t280;
t56 = pkin(2) * t57;
t53 = -t200 * t83 - t204 * t80;
t51 = -t200 * t80 + t204 * t83;
t48 = t202 * (-t109 * t201 + t110 * t205) + t206 * (t109 * t205 + t110 * t201);
t47 = t151 * qJ(5) - t292;
t43 = (-t136 - t151) * qJ(5) + t292;
t42 = t237 - t270;
t41 = t238 + t266;
t40 = t202 * (-t101 * t201 + t103 * t205) + t206 * (t101 * t205 + t103 * t201);
t39 = t202 * (-t100 * t201 + t102 * t205) + t206 * (t100 * t205 + t102 * t201);
t37 = t209 + (t113 + t151) * qJ(5) + (-t80 + t92) * pkin(4);
t30 = pkin(2) * t31;
t28 = t202 * (-t201 * t75 + t205 * t76) + t206 * (t201 * t76 + t205 * t75);
t27 = t202 * (-t201 * t73 + t205 * t74) + t206 * (t201 * t74 + t205 * t73);
t24 = pkin(3) * t25;
t23 = -t200 * t65 + t204 * t43 - t279;
t22 = (t84 + t93) * qJ(5) + t286;
t21 = qJ(5) * t251 - t200 * t37 - t280;
t20 = -pkin(4) * t106 + qJ(5) * t215 + t35;
t19 = pkin(3) * t90 + pkin(8) * t26;
t18 = t200 * t43 + t204 * t65 + t237;
t17 = pkin(4) * t47 + qJ(5) * t35;
t16 = qJ(5) * t258 + t204 * t37 + t238;
t15 = -t25 - t281;
t14 = t234 + t26;
t13 = t204 * t35 - t271;
t10 = t202 * (-t201 * t51 + t205 * t53) + t206 * (t201 * t53 + t205 * t51);
t8 = t205 * t26 - t269;
t7 = t201 * t26 + t265;
t6 = -t20 * t200 + t204 * t22 - t281;
t5 = t20 * t204 + t200 * t22 + t234;
t4 = -t12 * t201 + t13 * t205;
t3 = t12 * t205 + t13 * t201;
t2 = -pkin(8) * t12 - qJ(5) * t267 - t17 * t200;
t1 = pkin(3) * t47 + pkin(8) * t13 - qJ(5) * t271 + t17 * t204;
t9 = [0, 0, 0, 0, 0, qJDD(1), t231, t221, 0, 0, (t179 + t232) * t202, t178 * t206 + t181 * t202, t254 + t206 * (-t192 + t207), (t180 - t233) * t206, t202 * (t193 - t207) + t248, 0, t206 * t174 + pkin(1) * t181 + pkin(6) * (t206 * (-t193 - t207) - t254), -t202 * t174 - pkin(1) * t178 + pkin(6) * (-t248 - t202 * (-t192 - t207)), pkin(1) * (t192 + t193) + (t198 + t199) * t264 + t229, pkin(1) * t174 + pkin(6) * t229, t202 * (t138 * t205 - t173 * t261) + t206 * (t138 * t201 + t173 * t260), t202 * (-t122 * t205 - t126 * t201) + t206 * (-t122 * t201 + t126 * t205), t202 * (-t167 * t201 + t290) + t206 * (t167 * t205 + t291), t202 * (-t137 * t201 - t172 * t260) + t206 * (t137 * t205 - t172 * t261), t202 * (t166 * t205 - t256) + t206 * (t166 * t201 + t249), (t202 * (t172 * t205 + t173 * t201) + t206 * (t172 * t201 - t173 * t205)) * t197, t202 * (-pkin(7) * t117 - t257) + t206 * (-pkin(2) * t122 + pkin(7) * t118 + t250) - pkin(1) * t122 + pkin(6) * (-t117 * t202 + t118 * t206), t202 * (-pkin(7) * t128 - t250) + t206 * (-pkin(2) * t125 + pkin(7) * t129 - t257) - pkin(1) * t125 + pkin(6) * (-t128 * t202 + t129 * t206), t202 * (-pkin(7) * t94 - t66) + t206 * (-pkin(2) * t139 + pkin(7) * t95 + t67) - pkin(1) * t139 + pkin(6) * (-t202 * t94 + t206 * t95), -pkin(7) * t268 + t206 * (pkin(2) * t140 + pkin(7) * t67) + pkin(1) * t140 + pkin(6) * (t206 * t67 - t268), t28, t10, t39, t27, t40, t48, t202 * (-t201 * t41 + t205 * t59 - t283) + t206 * (t201 * t59 + t205 * t41 + t240) + t275, t202 * (-t201 * t42 + t205 * t60 - t282) + t206 * (t201 * t60 + t205 * t42 + t239) + t274, t202 * (-t14 * t201 + t15 * t205 - t284) + t206 * (t14 * t205 + t15 * t201 + t235) + t272, t202 * (-pkin(7) * t7 - pkin(8) * t265 - t19 * t201) + t206 * (pkin(2) * t90 + pkin(7) * t8 - pkin(8) * t269 + t19 * t205) + pkin(1) * t90 + pkin(6) * (-t202 * t7 + t206 * t8), t28, t10, t39, t27, t40, t48, t202 * (-t16 * t201 + t205 * t21 - t283) + t206 * (t16 * t205 + t201 * t21 + t240) + t275, t202 * (-t18 * t201 + t205 * t23 - t282) + t206 * (t18 * t205 + t201 * t23 + t239) + t274, t202 * (-t201 * t5 + t205 * t6 - t284) + t206 * (t201 * t6 + t205 * t5 + t235) + t272, t202 * (-pkin(7) * t3 - t1 * t201 + t2 * t205) + t206 * (pkin(2) * t47 + pkin(7) * t4 + t1 * t205 + t2 * t201) + pkin(1) * t47 + pkin(6) * (-t202 * t3 + t206 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185, -t193 + t192, t190, t185, t191, qJDD(2), -t164, -t165, 0, 0, -t161, t160, t127, t161, t214, t196, pkin(2) * t117 - t107, pkin(2) * t128 - t108, pkin(2) * t94, pkin(2) * t66, t121, t119, t84, -t121, t215, t188, t227 + t56, t62 + t241, t30 + t50, pkin(2) * t7 + t24, t121, t119, t84, -t121, t215, t188, t211 + t56, t62 + t216, t30 + t273, pkin(2) * t3 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, t160, t127, t161, t214, t196, -t107, -t108, 0, 0, t121, t119, t84, -t121, t215, t188, t227, t241, t50, t24, t121, t119, t84, -t121, t215, t188, t211, t216, t273, t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t119, t84, -t121, t215, t188, -t45, -t46, 0, 0, t121, t119, t84, -t121, t215, t188, t212, t217, -t77, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t85, t106, -t47;];
tauJ_reg = t9;
