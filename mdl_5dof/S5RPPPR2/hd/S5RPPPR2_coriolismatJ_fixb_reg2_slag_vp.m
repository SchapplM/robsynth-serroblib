% Calculate inertial parameters regressor of coriolis matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:50
% EndTime: 2022-01-23 08:59:54
% DurationCPUTime: 2.50s
% Computational Cost: add. (3350->237), mult. (7909->391), div. (0->0), fcn. (8655->8), ass. (0->209)
t189 = cos(pkin(9));
t190 = cos(pkin(8));
t192 = sin(qJ(5));
t279 = t192 * t190;
t187 = sin(pkin(8));
t193 = cos(qJ(5));
t287 = t187 * t193;
t156 = t189 * t287 - t279;
t309 = -t156 / 0.2e1;
t186 = sin(pkin(9));
t191 = cos(pkin(7));
t280 = t191 * t189;
t188 = sin(pkin(7));
t283 = t190 * t188;
t150 = t186 * t283 + t280;
t148 = t150 ^ 2;
t308 = t186 / 0.2e1;
t307 = -t189 / 0.2e1;
t306 = -t192 / 0.2e1;
t305 = -t193 / 0.2e1;
t276 = t193 * t191;
t286 = t188 * t186;
t154 = t190 * t280 + t286;
t298 = t154 * t192;
t113 = t187 * t276 - t298;
t278 = t192 * t191;
t230 = t187 * t278;
t297 = t154 * t193;
t114 = t230 + t297;
t281 = t191 * t186;
t285 = t189 * t188;
t149 = t190 * t281 - t285;
t152 = t189 * t283 - t281;
t214 = -t188 * qJ(3) - pkin(1);
t195 = -t190 * t214 + (t190 * pkin(2) + t187 * qJ(2) + pkin(3)) * t191;
t194 = t150 * pkin(4) - t152 * pkin(6) + t195;
t290 = t187 * t188;
t208 = -pkin(2) * t191 + t214;
t282 = t191 * qJ(2);
t130 = t187 * t208 + t190 * t282;
t115 = -t191 * qJ(4) + t130;
t143 = (pkin(3) * t187 - qJ(4) * t190 + qJ(2)) * t188;
t68 = t189 * t115 + t186 * t143;
t65 = pkin(6) * t290 + t68;
t31 = t192 * t65 - t193 * t194;
t32 = t192 * t194 + t193 * t65;
t67 = -t186 * t115 + t189 * t143;
t64 = -pkin(4) * t290 - t67;
t5 = -t31 * t113 + t32 * t114 + t64 * t149;
t304 = t5 * qJD(1);
t6 = t64 * t152 + (-t192 * t31 - t32 * t193) * t150;
t303 = t6 * qJD(1);
t288 = t187 * t192;
t231 = t189 * t288;
t277 = t193 * t190;
t155 = t231 + t277;
t131 = t155 * t188;
t132 = t156 * t188;
t233 = t187 * t286;
t7 = -t31 * t131 - t32 * t132 - t64 * t233;
t302 = t7 * qJD(1);
t110 = t192 * t152 - t188 * t287;
t12 = -t64 * t110 + t31 * t150;
t301 = t12 * qJD(1);
t112 = t152 * t193 + t188 * t288;
t13 = t64 * t112 - t32 * t150;
t300 = t13 * qJD(1);
t299 = t154 * t186;
t181 = t186 ^ 2;
t223 = t149 * t307;
t197 = (t114 * t193 / 0.2e1 + t113 * t306) * t186 + t223;
t206 = -t131 * t155 / 0.2e1 + t132 * t309;
t182 = t187 ^ 2;
t293 = t182 * t188;
t16 = t181 * t293 / 0.2e1 + t197 - t206;
t296 = t16 * qJD(1);
t289 = t187 * t191;
t18 = -t67 * t149 + t68 * t154 + t195 * t289;
t295 = t18 * qJD(1);
t183 = t188 ^ 2;
t294 = t182 * t183;
t292 = t186 * t150;
t291 = t186 * t187;
t232 = t187 * t285;
t19 = -t195 * t283 + t68 * t232 - t67 * t233;
t284 = t19 * qJD(1);
t216 = t291 / 0.2e1;
t196 = (t155 * t306 + t156 * t305) * t150 + t152 * t216;
t207 = t113 * t305 + t114 * t306;
t21 = t196 + t207;
t275 = t21 * qJD(1);
t26 = -t68 * t150 - t67 * t152;
t274 = t26 * qJD(1);
t33 = -t114 * t110 - t113 * t112;
t273 = t33 * qJD(1);
t222 = t152 * t307;
t198 = (-t193 ^ 2 / 0.2e1 - t192 ^ 2 / 0.2e1) * t292 + t222;
t205 = t131 * t305 - t132 * t306;
t35 = t198 + t205;
t272 = t35 * qJD(1);
t36 = t132 * t110 - t131 * t112;
t271 = t36 * qJD(1);
t218 = -t186 * t110 / 0.2e1;
t221 = t155 * t150 / 0.2e1;
t37 = t221 + t297 / 0.2e1 + (t218 + t278 / 0.2e1) * t187;
t270 = t37 * qJD(1);
t220 = t150 * t309;
t39 = t220 + t298 / 0.2e1 + (t112 * t308 - t276 / 0.2e1) * t187;
t269 = t39 * qJD(1);
t41 = (t193 * t110 - t112 * t192) * t150;
t268 = t41 * qJD(1);
t42 = t149 * t110 + t113 * t150;
t267 = t42 * qJD(1);
t43 = t149 * t112 - t114 * t150;
t266 = t43 * qJD(1);
t44 = t110 ^ 2 - t112 ^ 2;
t265 = t44 * qJD(1);
t45 = -t152 * t110 - t192 * t148;
t264 = t45 * qJD(1);
t200 = (t231 / 0.2e1 + t277 / 0.2e1) * t188;
t217 = -t292 / 0.2e1;
t202 = t112 * t307 + t193 * t217;
t46 = t200 - t202;
t263 = t46 * qJD(1);
t215 = t287 / 0.2e1;
t199 = (t189 * t215 - t279 / 0.2e1) * t188;
t201 = t110 * t307 + t192 * t217;
t47 = t199 + t201;
t262 = t47 * qJD(1);
t50 = -t110 * t233 + t131 * t150;
t261 = t50 * qJD(1);
t51 = -t112 * t233 + t132 * t150;
t260 = t51 * qJD(1);
t219 = -t293 / 0.2e1;
t163 = t181 * t219;
t184 = t190 ^ 2;
t174 = -t184 * t188 / 0.2e1;
t247 = t163 + t174;
t55 = -t299 / 0.2e1 + (t149 / 0.2e1 + t189 * t219) * t189 + t247;
t259 = t55 * qJD(1);
t57 = t149 * t152 - t154 * t150;
t258 = t57 * qJD(1);
t59 = t152 * t112 + t193 * t148;
t257 = t59 * qJD(1);
t129 = -t187 * t282 + t190 * t208;
t63 = t183 * qJ(2) + (-t129 * t187 + t130 * t190) * t191;
t256 = t63 * qJD(1);
t66 = (t129 * t190 + t130 * t187) * t188;
t255 = t66 * qJD(1);
t203 = t150 * t307 + t152 * t308;
t70 = (-t191 / 0.2e1 + t203) * t187;
t254 = t70 * qJD(1);
t71 = -t150 * t232 + t152 * t233;
t253 = t71 * qJD(1);
t204 = t222 + t217;
t73 = -t283 / 0.2e1 + t204;
t252 = t73 * qJD(1);
t79 = (-t149 * t188 + t150 * t191) * t187;
t251 = t79 * qJD(1);
t80 = (t152 * t191 - t154 * t188) * t187;
t250 = t80 * qJD(1);
t84 = t192 * t150;
t249 = t84 * qJD(1);
t97 = t152 ^ 2 + t148;
t248 = t97 * qJD(1);
t171 = t191 ^ 2 + t183;
t246 = qJD(1) * t150;
t245 = qJD(1) * t188;
t244 = qJD(5) * t192;
t243 = qJD(5) * t193;
t107 = t150 * t283 + t186 * t294;
t242 = t107 * qJD(1);
t108 = t152 * t283 + t189 * t294;
t241 = t108 * qJD(1);
t240 = t112 * qJD(5);
t146 = (0.1e1 / 0.2e1 + t182 / 0.2e1 + t184 / 0.2e1) * t188;
t239 = t146 * qJD(1);
t158 = (t182 + t184) * t183;
t238 = t158 * qJD(1);
t159 = t171 * t187;
t237 = t159 * qJD(1);
t160 = t171 * t190;
t236 = t160 * qJD(1);
t169 = t171 * qJ(2);
t235 = t169 * qJD(1);
t234 = t171 * qJD(1);
t229 = t187 * t245;
t228 = t191 * t245;
t227 = qJD(3) * t188 * t191;
t226 = qJD(4) * t290;
t225 = t112 * t110 * qJD(1);
t224 = t110 * t240;
t213 = -qJD(5) - t246;
t212 = t152 * t229;
t211 = t150 * t229;
t210 = t187 * t228;
t209 = t190 * t228;
t147 = t219 + t174 + t188 / 0.2e1;
t72 = t283 / 0.2e1 + t204;
t69 = t289 / 0.2e1 + t203 * t187;
t54 = t189 ^ 2 * t219 + t299 / 0.2e1 + t223 + t247;
t49 = t200 + t202;
t48 = t199 - t201;
t40 = t220 + t112 * t216 - t298 / 0.2e1 + t191 * t215;
t38 = t221 + t187 * t218 - t297 / 0.2e1 - t230 / 0.2e1;
t34 = t198 - t205;
t20 = t196 - t207;
t17 = t163 + t197 + t206;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171 * qJD(2), t169 * qJD(2), 0, 0, 0, 0, 0, 0, t159 * qJD(2) + t190 * t227, t160 * qJD(2) - t187 * t227, t158 * qJD(3), t63 * qJD(2) - t66 * qJD(3), 0, 0, 0, 0, 0, 0, t79 * qJD(2) + t107 * qJD(3) - t152 * t226, t80 * qJD(2) + t108 * qJD(3) + t150 * t226, t57 * qJD(2) - t71 * qJD(3) + t97 * qJD(4), t18 * qJD(2) - t19 * qJD(3) + t26 * qJD(4), -t224, t44 * qJD(5), -t110 * t150 * qJD(5), t224, -t150 * t240, 0, t42 * qJD(2) + t50 * qJD(3) - t45 * qJD(4) + t13 * qJD(5), t43 * qJD(2) + t51 * qJD(3) + t59 * qJD(4) + t12 * qJD(5), t33 * qJD(2) + t36 * qJD(3) + t41 * qJD(4), t5 * qJD(2) + t7 * qJD(3) + t6 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, t235, 0, 0, 0, 0, 0, 0, t237, t236, 0, t147 * qJD(3) + t256, 0, 0, 0, 0, 0, 0, t251, t250, t258, t295 + t54 * qJD(3) + t69 * qJD(4) + (t149 * t186 + t154 * t189 - t190 * t191) * qJD(2) * t187, 0, 0, 0, 0, 0, 0, t40 * qJD(5) + t267, t38 * qJD(5) + t266, t273, t304 + (-t113 * t155 + t114 * t156 + t149 * t291) * qJD(2) + t17 * qJD(3) + t20 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, -t210, t238, t147 * qJD(2) - t255, 0, 0, 0, 0, 0, 0, t242, t241, -t253, t54 * qJD(2) + t72 * qJD(4) - t284, 0, 0, 0, 0, 0, 0, t49 * qJD(5) + t261, t48 * qJD(5) + t260, t271, t302 + t17 * qJD(2) + t34 * qJD(4) + (-t131 * t192 - t132 * t193 + t232) * qJD(3) * t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212, t211, t248, t69 * qJD(2) + t72 * qJD(3) + t274, 0, 0, 0, 0, 0, 0, -t264, t257, t268, t20 * qJD(2) + t34 * qJD(3) + t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225, t265, t213 * t110, t225, t213 * t112, 0, t40 * qJD(2) + t49 * qJD(3) - t32 * qJD(5) + t300, t38 * qJD(2) + t48 * qJD(3) + t31 * qJD(5) + t301, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, -t235, 0, 0, 0, 0, 0, 0, -t237, -t236, 0, -t146 * qJD(3) - t256, 0, 0, 0, 0, 0, 0, -t251, -t250, -t258, t55 * qJD(3) + t70 * qJD(4) - t295, 0, 0, 0, 0, 0, 0, t39 * qJD(5) - t267, t37 * qJD(5) - t266, -t273, -t16 * qJD(3) + t21 * qJD(4) - t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156 * qJD(5) + t269, t155 * qJD(5) + t270, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t210, -t238, t146 * qJD(2) + t255, 0, 0, 0, 0, 0, 0, -t242, -t241, t253, -t55 * qJD(2) + t73 * qJD(4) + t284, 0, 0, 0, 0, 0, 0, -t46 * qJD(5) - t261, -t47 * qJD(5) - t260, -t271, t16 * qJD(2) + t35 * qJD(4) - t302; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t239, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t186 * t243 - t263, t186 * t244 - t262, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, -t211, -t248, -t70 * qJD(2) - t73 * qJD(3) - t274, 0, 0, 0, 0, 0, 0, -t84 * qJD(5) + t264, -t150 * t243 - t257, -t268, -t21 * qJD(2) - t35 * qJD(3) - t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t252, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t272; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244 - t249, t213 * t193, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, -t265, t110 * t246, -t225, t112 * t246, 0, -t39 * qJD(2) + t46 * qJD(3) + t84 * qJD(4) - t300, t193 * t150 * qJD(4) - t37 * qJD(2) + t47 * qJD(3) - t301, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t269, -t270, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, t262, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, t193 * t246, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
