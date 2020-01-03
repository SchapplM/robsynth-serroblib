% Calculate inertial parameters regressor of coriolis matrix for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR7_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:56
% EndTime: 2019-12-31 18:20:00
% DurationCPUTime: 2.64s
% Computational Cost: add. (3782->208), mult. (7111->300), div. (0->0), fcn. (7525->8), ass. (0->182)
t173 = sin(qJ(5));
t174 = sin(qJ(3));
t165 = sin(pkin(8)) * pkin(1) + pkin(6);
t289 = qJ(4) + t165;
t148 = t289 * t174;
t172 = sin(pkin(9));
t176 = cos(qJ(3));
t238 = t289 * t176;
t264 = cos(pkin(9));
t84 = t264 * t148 + t172 * t238;
t291 = t84 * t173;
t175 = cos(qJ(5));
t290 = t84 * t175;
t286 = -t172 * t148 + t264 * t238;
t269 = t286 * t173;
t268 = t286 * t175;
t202 = t264 * t174;
t256 = t172 * t176;
t154 = t202 + t256;
t254 = t173 * t175;
t201 = 0.2e1 * t154 * t254;
t147 = t154 * qJD(3);
t140 = t175 * t147;
t152 = t172 * t174 - t264 * t176;
t98 = t173 * t152;
t91 = t98 * qJD(5);
t288 = -t140 + t91;
t150 = t152 ^ 2;
t151 = t154 ^ 2;
t287 = -t151 - t150;
t217 = t151 - t150;
t285 = t286 / 0.2e1;
t284 = t84 / 0.2e1;
t171 = t175 ^ 2;
t283 = -t171 / 0.2e1;
t282 = t174 * pkin(3);
t280 = qJD(3) * pkin(3);
t103 = t154 * pkin(4) + t152 * pkin(7) + t282;
t255 = t173 * t103;
t53 = t255 - t290;
t272 = t53 * t173;
t253 = t175 * t103;
t52 = t253 + t291;
t273 = t52 * t175;
t167 = -cos(pkin(8)) * pkin(1) - pkin(2);
t156 = -t176 * pkin(3) + t167;
t177 = t152 * pkin(4) - t154 * pkin(7) + t156;
t48 = -t175 * t177 + t269;
t49 = t173 * t177 + t268;
t4 = (t272 + t273) * t154 + (-t49 * t173 + t175 * t48) * t152;
t277 = t4 * qJD(1);
t276 = t48 * t173;
t275 = t49 * t175;
t274 = t52 * t173;
t271 = t53 * t175;
t6 = (-t48 + t269) * t154 + (t52 - t291) * t152;
t270 = t6 * qJD(1);
t267 = t84 * t154;
t170 = t173 ^ 2;
t204 = t283 - t170 / 0.2e1;
t164 = t172 * pkin(3) + pkin(7);
t258 = t164 * t152;
t166 = -t264 * pkin(3) - pkin(4);
t259 = t154 * t166;
t178 = t204 * t258 + t259 / 0.2e1;
t188 = -t273 / 0.2e1 - t272 / 0.2e1;
t11 = t178 + t188;
t263 = t11 * qJD(1);
t216 = -t84 / 0.2e1 + t284;
t12 = t216 * t154 + (-t286 / 0.2e1 + t285) * t152;
t262 = t12 * qJD(1);
t261 = t152 * t154;
t260 = t152 * t172;
t100 = t173 * t154;
t28 = -t84 * t100 + t48 * t152;
t252 = t28 * qJD(1);
t29 = -t49 * t152 + t154 * t290;
t251 = t29 * qJD(1);
t30 = -t152 * t286 + t267;
t250 = t30 * qJD(1);
t40 = (0.1e1 / 0.2e1 + t204) * t261;
t248 = t40 * qJD(1);
t65 = t217 * t173;
t246 = t65 * qJD(1);
t66 = t287 * t173;
t245 = t66 * qJD(1);
t67 = t217 * t175;
t244 = t67 * qJD(1);
t203 = t264 * t154;
t180 = -t260 / 0.2e1 - t203 / 0.2e1;
t73 = (-t174 / 0.2e1 + t180) * pkin(3);
t243 = t73 * qJD(1);
t74 = t152 * t282 + t156 * t154;
t242 = t74 * qJD(1);
t75 = -t156 * t152 + t154 * t282;
t241 = t75 * qJD(1);
t240 = t217 * qJD(1);
t96 = (t170 / 0.2e1 + t283) * t154;
t239 = t96 * qJD(5);
t90 = t98 * qJD(1);
t237 = -t170 - t171;
t160 = t171 - t170;
t236 = qJD(1) * t176;
t235 = qJD(3) * t175;
t234 = qJD(5) * t173;
t233 = qJD(5) * t175;
t232 = t100 * qJD(1);
t102 = t175 * t152;
t231 = t102 * qJD(1);
t230 = t102 * qJD(3);
t142 = t170 * t152;
t143 = t171 * t152;
t104 = t142 + t143;
t229 = t104 * qJD(1);
t228 = t104 * qJD(3);
t106 = t287 * t175;
t227 = t106 * qJD(1);
t226 = t287 * qJD(1);
t149 = t202 / 0.2e1 + t256 / 0.2e1;
t225 = t149 * qJD(1);
t224 = t152 * qJD(1);
t146 = t152 * qJD(3);
t223 = t152 * qJD(4);
t222 = t154 * qJD(1);
t221 = t154 * qJD(4);
t161 = -t174 ^ 2 + t176 ^ 2;
t220 = t161 * qJD(1);
t219 = t174 * qJD(3);
t218 = t176 * qJD(3);
t214 = t154 * t234;
t213 = t154 * t233;
t212 = t152 * t222;
t211 = t152 * t147;
t210 = t167 * t174 * qJD(1);
t209 = t167 * t236;
t208 = t173 * t233;
t207 = t173 * t235;
t206 = t174 * t218;
t205 = t175 * t222;
t200 = -qJD(5) - t224;
t199 = t151 * t208;
t198 = qJD(3) * t201;
t2 = (t271 / 0.2e1 - t274 / 0.2e1 + t284) * t154 + (-t275 / 0.2e1 - t276 / 0.2e1 + t285) * t152;
t3 = t84 * t286 - t48 * t52 + t49 * t53;
t197 = t3 * qJD(1) + t2 * qJD(2);
t196 = t271 - t274;
t46 = (0.1e1 + t237) * t261;
t195 = -t2 * qJD(1) - t46 * qJD(2);
t7 = (-t49 + t268) * t154 + (-t53 - t290) * t152;
t194 = t7 * qJD(1);
t9 = t267 + (-t275 - t276) * t152;
t193 = t9 * qJD(1) + t40 * qJD(2);
t26 = t156 * t282;
t192 = t26 * qJD(1) + t12 * qJD(2);
t190 = -t152 * t166 - t154 * t164;
t189 = t200 * t175;
t187 = t258 / 0.2e1 - t259 / 0.2e1;
t179 = t103 / 0.2e1 + t187;
t24 = t216 * t173 - t179 * t175;
t186 = -t166 * t173 * qJD(3) - t24 * qJD(1);
t22 = t179 * t173 + t216 * t175;
t185 = -t22 * qJD(1) - t166 * t235;
t69 = -t96 * qJD(1) + t207;
t184 = t154 * t189;
t183 = t149 * qJD(5) + t212;
t55 = t151 * qJD(1) * t254 + t96 * qJD(3);
t105 = t160 * t151;
t182 = t105 * qJD(1) + t198;
t181 = qJD(1) * t201 - t160 * qJD(3);
t162 = t174 * t236;
t141 = t149 * qJD(3);
t139 = t173 * t147;
t95 = t102 * qJD(5);
t89 = t98 * qJD(3);
t72 = t282 / 0.2e1 + t180 * pkin(3);
t71 = -t234 - t90;
t25 = t291 + t253 / 0.2e1 - t187 * t175;
t23 = t290 - t255 / 0.2e1 + t187 * t173;
t10 = t178 - t188;
t8 = t12 * qJD(3);
t1 = t2 * qJD(3) + t40 * qJD(4);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, t161 * qJD(3), 0, -t206, 0, 0, t167 * t219, t167 * t218, 0, 0, -t211, -t217 * qJD(3), 0, t211, 0, 0, t74 * qJD(3), t75 * qJD(3), -qJD(4) * t287, t26 * qJD(3) + t30 * qJD(4), -t171 * t211 - t199, -t105 * qJD(5) + t152 * t198, t67 * qJD(3) - t152 * t214, -t170 * t211 + t199, -t65 * qJD(3) - t152 * t213, t211, t6 * qJD(3) - t66 * qJD(4) + t29 * qJD(5), t7 * qJD(3) - t106 * qJD(4) + t28 * qJD(5), -t4 * qJD(3), t3 * qJD(3) + t9 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t220, t218, -t162, -t219, 0, -t165 * t218 + t210, t165 * t219 + t209, 0, 0, -t212, -t240, -t146, t212, -t147, 0, -qJD(3) * t286 + t242, qJD(3) * t84 + t241, (t152 * t264 - t154 * t172) * t280, (-t172 * t84 - t264 * t286) * t280 + t72 * qJD(4) + t192, -t239 + (-t171 * t222 - t207) * t152, (t142 - t143) * qJD(3) + (-qJD(5) + t224) * t201, t139 + t244, t239 + (-t170 * t222 + t207) * t152, t140 - t246, t183, t270 + (t173 * t190 - t268) * qJD(3) + t25 * qJD(5), (t175 * t190 + t269) * qJD(3) + t23 * qJD(5) + t194, qJD(3) * t196 - t277, (t164 * t196 + t166 * t286) * qJD(3) + t10 * qJD(4) + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t226, t72 * qJD(3) + t250, 0, 0, 0, 0, 0, 0, -t245, -t227, 0, t10 * qJD(3) + t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t182, t200 * t100, t55, t184, t141, t25 * qJD(3) - t49 * qJD(5) + t251, t23 * qJD(3) + t48 * qJD(5) + t252, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, -t218, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t146, 0, t262 + (-t203 - t260) * t280, 0, 0, 0, 0, 0, 0, t288, t139 + t95, -t228, (t237 * t258 + t259) * qJD(3) - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 - t213, t214 + t230, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, -t220, 0, t162, 0, 0, -t210, -t209, 0, 0, t212, t240, 0, -t212, 0, 0, -t221 - t242, t223 - t241, 0, t73 * qJD(4) - t192, t171 * t212 - t239, 0.2e1 * t173 * t184, t95 - t244, t170 * t212 + t239, -t91 + t246, -t183, t24 * qJD(5) - t175 * t221 - t270, t100 * qJD(4) + t22 * qJD(5) - t194, -t104 * qJD(4) + t277, t11 * qJD(4) - t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t160 * qJD(5), 0, -t208, 0, 0, t166 * t234, t166 * t233, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222, t224, 0, t243, 0, 0, 0, 0, 0, 0, -t205, t232, -t229, t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t181, t231 + t233, -t69, t71, -t225, -t164 * t233 - t186, t164 * t234 - t185, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, -t146, t226, -t73 * qJD(3) - t250, 0, 0, 0, 0, 0, 0, t245 - t288, -t100 * qJD(3) - t152 * t233 + t227, t228, -t11 * qJD(3) - t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, -t224, 0, -t243, 0, 0, 0, 0, 0, 0, t205, -t232, t229, -t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t189, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t182, t173 * t212 - t230, -t55, t152 * t205 + t89, t141, -t24 * qJD(3) + t98 * qJD(4) - t251, -t22 * qJD(3) + t175 * t223 - t252, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t181, -t231, t69, t90, t225, t186, t185, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t175 * t224, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
