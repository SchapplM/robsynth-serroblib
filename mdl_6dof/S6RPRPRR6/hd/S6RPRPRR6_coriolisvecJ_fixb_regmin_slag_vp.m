% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:53
% EndTime: 2019-03-09 03:53:05
% DurationCPUTime: 4.80s
% Computational Cost: add. (6783->386), mult. (17950->533), div. (0->0), fcn. (14732->10), ass. (0->187)
t208 = cos(pkin(10));
t285 = cos(qJ(3));
t245 = t285 * t208;
t197 = qJD(1) * t245;
t206 = sin(pkin(10));
t211 = sin(qJ(3));
t264 = t211 * t206;
t244 = qJD(1) * t264;
t166 = -t197 + t244;
t160 = qJD(5) + t166;
t154 = qJD(6) + t160;
t209 = sin(qJ(6));
t212 = cos(qJ(6));
t183 = t285 * t206 + t211 * t208;
t168 = t183 * qJD(1);
t205 = sin(pkin(11));
t207 = cos(pkin(11));
t143 = -t207 * qJD(3) + t205 * t168;
t145 = t205 * qJD(3) + t207 * t168;
t210 = sin(qJ(5));
t213 = cos(qJ(5));
t290 = t210 * t143 - t213 * t145;
t94 = t213 * t143 + t210 * t145;
t41 = -t209 * t290 + t212 * t94;
t304 = t41 * t154;
t303 = t94 * t160;
t229 = t209 * t94 + t212 * t290;
t302 = t154 * t229;
t262 = t213 * t207;
t180 = t210 * t205 - t262;
t261 = t160 * t180;
t182 = t213 * t205 + t210 * t207;
t170 = t182 * qJD(5);
t260 = t182 * t166 + t170;
t271 = t166 * t205;
t129 = t168 * pkin(3) + t166 * qJ(4);
t283 = pkin(7) + qJ(2);
t188 = t283 * t206;
t184 = qJD(1) * t188;
t190 = t283 * t208;
t185 = qJD(1) * t190;
t222 = t285 * t184 + t211 * t185;
t78 = t205 * t129 - t207 * t222;
t62 = pkin(8) * t271 + t78;
t301 = -qJD(4) * t207 + t62;
t300 = t160 * t290;
t133 = -t209 * t180 + t212 * t182;
t278 = qJD(6) * t133 - t261 * t209 + t260 * t212;
t297 = t229 * t41;
t296 = t183 * qJD(2);
t295 = t229 ^ 2 - t41 ^ 2;
t201 = -t208 * pkin(2) - pkin(1);
t186 = qJD(1) * t201 + qJD(2);
t110 = t166 * pkin(3) - t168 * qJ(4) + t186;
t135 = -t211 * t184 + t285 * t185;
t128 = qJD(3) * qJ(4) + t135;
t67 = t207 * t110 - t205 * t128;
t34 = t166 * pkin(4) - t145 * pkin(8) + t67;
t68 = t205 * t110 + t207 * t128;
t47 = -t143 * pkin(8) + t68;
t19 = t210 * t34 + t213 * t47;
t13 = -t94 * pkin(9) + t19;
t252 = qJD(6) * t209;
t11 = t13 * t252;
t126 = -qJD(3) * pkin(3) + qJD(4) + t222;
t84 = t143 * pkin(4) + t126;
t39 = t94 * pkin(5) + t84;
t294 = t39 * t41 + t11;
t251 = qJD(6) * t212;
t193 = qJD(3) * t197;
t155 = -qJD(3) * t244 + t193;
t253 = qJD(5) * t213;
t266 = t205 * t155;
t52 = -t143 * t253 + t155 * t262 + (-qJD(5) * t145 - t266) * t210;
t53 = -qJD(5) * t290 + t155 * t182;
t8 = -t209 * t53 + t212 * t52 - t94 * t251 + t252 * t290;
t293 = t8 + t304;
t172 = t183 * qJD(3);
t156 = qJD(1) * t172;
t265 = t207 * t155;
t93 = t156 * pkin(3) - t155 * qJ(4) - t168 * qJD(4);
t221 = t245 - t264;
t219 = t221 * qJD(2);
t99 = qJD(1) * t219 + (qJD(4) - t222) * qJD(3);
t37 = -t205 * t99 + t207 * t93;
t25 = t156 * pkin(4) - pkin(8) * t265 + t37;
t38 = t205 * t93 + t207 * t99;
t28 = -pkin(8) * t266 + t38;
t239 = -t210 * t28 + t213 * t25;
t216 = -qJD(5) * t19 + t239;
t2 = t156 * pkin(5) - t52 * pkin(9) + t216;
t254 = qJD(5) * t210;
t224 = t210 * t25 + t213 * t28 + t34 * t253 - t254 * t47;
t3 = -t53 * pkin(9) + t224;
t249 = t212 * t2 - t209 * t3;
t18 = -t210 * t47 + t213 * t34;
t12 = pkin(9) * t290 + t18;
t10 = t160 * pkin(5) + t12;
t275 = t212 * t13;
t5 = t209 * t10 + t275;
t292 = -qJD(6) * t5 + t39 * t229 + t249;
t215 = qJD(6) * t229 - t209 * t52 - t212 * t53;
t291 = t215 - t302;
t282 = pkin(8) + qJ(4);
t187 = t282 * t205;
t189 = t282 * t207;
t258 = -t210 * t187 + t213 * t189;
t289 = -t285 * t188 - t211 * t190;
t132 = t212 * t180 + t209 * t182;
t279 = -qJD(6) * t132 - t260 * t209 - t261 * t212;
t288 = -t133 * t156 - t279 * t154;
t287 = -t182 * t156 + t160 * t261;
t225 = qJD(4) * t205 + qJD(5) * t189;
t284 = pkin(8) * t207;
t77 = t207 * t129 + t205 * t222;
t48 = t168 * pkin(4) + t166 * t284 + t77;
t286 = t187 * t253 + t301 * t213 + (t225 + t48) * t210;
t161 = t166 ^ 2;
t267 = t183 * t207;
t130 = -pkin(3) * t221 - t183 * qJ(4) + t201;
t142 = -t211 * t188 + t285 * t190;
t81 = t207 * t130 - t205 * t142;
t59 = -pkin(4) * t221 - pkin(8) * t267 + t81;
t268 = t183 * t205;
t82 = t205 * t130 + t207 * t142;
t70 = -pkin(8) * t268 + t82;
t280 = t210 * t59 + t213 * t70;
t277 = t168 * t41;
t276 = t168 * t94;
t274 = t229 * t168;
t273 = t290 * t168;
t171 = t221 * qJD(3);
t270 = t171 * t205;
t105 = t172 * pkin(3) - t171 * qJ(4) - t183 * qJD(4);
t111 = qJD(3) * t289 + t219;
t58 = t205 * t105 + t207 * t111;
t257 = t206 ^ 2 + t208 ^ 2;
t256 = qJD(3) * t211;
t250 = qJD(1) * qJD(2);
t103 = -pkin(4) * t271 + t135;
t200 = -t207 * pkin(4) - pkin(3);
t243 = t260 * pkin(5) - t103;
t242 = qJD(3) * t285;
t241 = qJD(6) * t10 + t3;
t57 = t207 * t105 - t205 * t111;
t31 = t172 * pkin(4) - t171 * t284 + t57;
t36 = -pkin(8) * t270 + t58;
t238 = -t210 * t36 + t213 * t31;
t237 = -t210 * t70 + t213 * t59;
t236 = t257 * qJD(1) ^ 2;
t235 = -t213 * t187 - t210 * t189;
t102 = t296 * qJD(1) - t184 * t256 + t185 * t242;
t112 = -t188 * t256 + t190 * t242 + t296;
t234 = -t132 * t156 - t278 * t154;
t117 = -t182 * pkin(9) + t235;
t233 = t260 * pkin(9) - qJD(6) * t117 + t286;
t118 = -t180 * pkin(9) + t258;
t46 = t213 * t48;
t232 = t168 * pkin(5) - t261 * pkin(9) + t182 * qJD(4) + t258 * qJD(5) + qJD(6) * t118 - t210 * t62 + t46;
t231 = -t180 * t156 - t260 * t160;
t230 = -t205 * t67 + t207 * t68;
t76 = pkin(4) * t266 + t102;
t83 = pkin(4) * t270 + t112;
t122 = t182 * t183;
t123 = t180 * t183;
t73 = t212 * t122 - t209 * t123;
t74 = -t209 * t122 - t212 * t123;
t226 = 0.2e1 * t257 * t250;
t113 = pkin(4) * t268 - t289;
t223 = t210 * t31 + t213 * t36 + t59 * t253 - t254 * t70;
t220 = t102 * t183 + t126 * t171 - t155 * t289;
t218 = -pkin(3) * t155 - qJ(4) * t156 + (-qJD(4) + t126) * t166;
t149 = t180 * pkin(5) + t200;
t136 = t156 * t221;
t75 = t122 * pkin(5) + t113;
t72 = t171 * t182 + t253 * t267 - t254 * t268;
t71 = -t170 * t183 - t171 * t180;
t26 = t72 * pkin(5) + t83;
t22 = t53 * pkin(5) + t76;
t21 = -t122 * pkin(9) + t280;
t20 = -pkin(5) * t221 + t123 * pkin(9) + t237;
t16 = qJD(6) * t74 + t209 * t71 + t212 * t72;
t15 = -qJD(6) * t73 - t209 * t72 + t212 * t71;
t7 = -t72 * pkin(9) + t223;
t6 = t172 * pkin(5) - t71 * pkin(9) - qJD(5) * t280 + t238;
t4 = t212 * t10 - t209 * t13;
t1 = [0, 0, 0, 0, 0, t226, qJ(2) * t226, t155 * t183 + t168 * t171, t155 * t221 - t183 * t156 - t171 * t166 - t168 * t172, t171 * qJD(3), -t172 * qJD(3), 0, -t112 * qJD(3) + t201 * t156 + t186 * t172, -t111 * qJD(3) + t201 * t155 + t186 * t171, t112 * t143 + t81 * t156 + t57 * t166 + t67 * t172 + t205 * t220 - t221 * t37, t112 * t145 - t82 * t156 - t58 * t166 - t68 * t172 + t207 * t220 + t221 * t38, -t58 * t143 - t57 * t145 + (-t155 * t81 - t171 * t67 - t183 * t37) * t207 + (-t155 * t82 - t171 * t68 - t183 * t38) * t205, -t102 * t289 + t126 * t112 + t37 * t81 + t38 * t82 + t67 * t57 + t68 * t58, -t52 * t123 - t290 * t71, -t52 * t122 + t123 * t53 + t290 * t72 - t71 * t94, -t123 * t156 + t71 * t160 - t172 * t290 - t221 * t52, -t122 * t156 - t72 * t160 - t94 * t172 + t221 * t53, t160 * t172 - t136, t238 * t160 + t237 * t156 - t239 * t221 + t18 * t172 + t83 * t94 + t113 * t53 + t76 * t122 + t84 * t72 + (-t160 * t280 + t19 * t221) * qJD(5), t113 * t52 - t76 * t123 - t280 * t156 - t223 * t160 - t19 * t172 + t221 * t224 - t290 * t83 + t84 * t71, -t15 * t229 + t8 * t74, -t15 * t41 + t16 * t229 + t215 * t74 - t8 * t73, t15 * t154 + t74 * t156 - t172 * t229 - t221 * t8, -t16 * t154 - t73 * t156 - t41 * t172 - t215 * t221, t154 * t172 - t136 (-t209 * t7 + t212 * t6) * t154 + (t212 * t20 - t209 * t21) * t156 - t249 * t221 + t4 * t172 + t26 * t41 - t75 * t215 + t22 * t73 + t39 * t16 + ((-t209 * t20 - t212 * t21) * t154 + t5 * t221) * qJD(6), -t11 * t221 + t39 * t15 - t5 * t172 + t22 * t74 - t26 * t229 + t75 * t8 + (-(-qJD(6) * t21 + t6) * t154 - t20 * t156 + t2 * t221) * t209 + (-(qJD(6) * t20 + t7) * t154 - t21 * t156 + t241 * t221) * t212; 0, 0, 0, 0, 0, -t236, -qJ(2) * t236, 0, 0, 0, 0, 0, 0.2e1 * t168 * qJD(3), t193 + (-t166 - t244) * qJD(3), -t168 * t143 + t207 * t156 - t205 * t161, -t168 * t145 - t205 * t156 - t207 * t161 (-t143 * t207 + t145 * t205) * t166 + (-t205 ^ 2 - t207 ^ 2) * t155, -t126 * t168 + t166 * t230 + t38 * t205 + t37 * t207, 0, 0, 0, 0, 0, t231 - t276, t273 + t287, 0, 0, 0, 0, 0, t234 - t277, t274 + t288; 0, 0, 0, 0, 0, 0, 0, t168 * t166, t168 ^ 2 - t161, t193 + (t166 - t244) * qJD(3), 0, 0, t135 * qJD(3) - t186 * t168 - t102, t186 * t166 - t221 * t250, -t102 * t207 - t135 * t143 - t77 * t166 - t67 * t168 + t205 * t218, t102 * t205 - t135 * t145 + t78 * t166 + t68 * t168 + t207 * t218, t78 * t143 + t77 * t145 + (-qJD(4) * t143 - t166 * t67 + t38) * t207 + (qJD(4) * t145 - t166 * t68 - t37) * t205, -t102 * pkin(3) - t126 * t135 - t67 * t77 - t68 * t78 + t230 * qJD(4) + (-t37 * t205 + t38 * t207) * qJ(4), t52 * t182 + t261 * t290, -t52 * t180 - t182 * t53 + t260 * t290 + t261 * t94, t273 - t287, t231 + t276, -t160 * t168, t235 * t156 + t200 * t53 + t76 * t180 - t18 * t168 - t103 * t94 + t260 * t84 + (-t46 - t225 * t213 + (qJD(5) * t187 + t301) * t210) * t160, t103 * t290 - t258 * t156 + t160 * t286 + t19 * t168 + t76 * t182 + t200 * t52 - t261 * t84, t8 * t133 - t229 * t279, -t8 * t132 + t133 * t215 + t229 * t278 - t279 * t41, t274 - t288, t234 + t277, -t154 * t168 (t212 * t117 - t209 * t118) * t156 - t149 * t215 + t22 * t132 - t4 * t168 + t243 * t41 + t278 * t39 + (t209 * t233 - t212 * t232) * t154 -(t209 * t117 + t212 * t118) * t156 + t149 * t8 + t22 * t133 + t5 * t168 - t243 * t229 + t279 * t39 + (t209 * t232 + t212 * t233) * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145 * t166 + t266, -t143 * t166 + t265, -t143 ^ 2 - t145 ^ 2, t143 * t68 + t145 * t67 + t102, 0, 0, 0, 0, 0, t53 - t300, t52 - t303, 0, 0, 0, 0, 0, -t215 - t302, t8 - t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290 * t94, t290 ^ 2 - t94 ^ 2, t52 + t303, -t53 - t300, t156, t19 * t160 + t290 * t84 + t216, t18 * t160 + t84 * t94 - t224, -t297, t295, t293, t291, t156 -(-t209 * t12 - t275) * t154 + (-t154 * t252 + t212 * t156 + t290 * t41) * pkin(5) + t292 (-t13 * t154 - t2) * t209 + (t12 * t154 - t241) * t212 + (-t154 * t251 - t209 * t156 - t229 * t290) * pkin(5) + t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t297, t295, t293, t291, t156, t5 * t154 + t292, t4 * t154 - t209 * t2 - t212 * t241 + t294;];
tauc_reg  = t1;
