% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_6_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:22
% EndTime: 2018-12-10 18:38:23
% DurationCPUTime: 1.69s
% Computational Cost: add. (20393->112), mult. (20452->213), div. (115->9), fcn. (20918->33), ass. (0->118)
t240 = sin(qJ(2));
t241 = sin(qJ(1));
t246 = cos(qJ(1));
t271 = pkin(6) + qJ(2);
t258 = cos(t271) / 0.2e1;
t272 = pkin(6) - qJ(2);
t263 = cos(t272);
t254 = t263 / 0.2e1 + t258;
t203 = t241 * t240 - t246 * t254;
t256 = sin(t271) / 0.2e1;
t261 = sin(t272);
t214 = t256 - t261 / 0.2e1;
t245 = cos(qJ(2));
t204 = t246 * t214 + t241 * t245;
t231 = pkin(7) - pkin(14);
t219 = cos(t231) / 0.2e1;
t230 = pkin(7) + pkin(14);
t228 = cos(t230);
t209 = t219 + t228 / 0.2e1;
t232 = sin(pkin(14));
t234 = sin(pkin(6));
t218 = sin(t230) / 0.2e1;
t227 = sin(t231);
t268 = t218 + t227 / 0.2e1;
t259 = t234 * t268;
t191 = t203 * t209 + t204 * t232 + t246 * t259;
t208 = t218 - t227 / 0.2e1;
t235 = cos(pkin(14));
t210 = t219 - t228 / 0.2e1;
t273 = t234 * t210;
t192 = t203 * t208 - t204 * t235 + t246 * t273;
t270 = pkin(8) - qJ(4);
t260 = sin(t270);
t269 = pkin(8) + qJ(4);
t286 = sin(t269) / 0.2e1;
t212 = t286 - t260 / 0.2e1;
t257 = cos(t269) / 0.2e1;
t262 = cos(t270);
t215 = t257 - t262 / 0.2e1;
t244 = cos(qJ(4));
t233 = sin(pkin(7));
t285 = cos(pkin(7));
t264 = t285 * t234;
t253 = t203 * t233 - t246 * t264;
t174 = t191 * t212 + t192 * t244 + t215 * t253;
t238 = sin(qJ(5));
t243 = cos(qJ(5));
t283 = sin(pkin(8));
t284 = cos(pkin(8));
t249 = t191 * t283 + t253 * t284;
t162 = t174 * t243 - t238 * t249;
t288 = t174 * t238 + t243 * t249;
t211 = t286 + t260 / 0.2e1;
t216 = t262 / 0.2e1 + t257;
t239 = sin(qJ(4));
t287 = -t191 * t216 + t192 * t239 + t211 * t253;
t213 = t256 + t261 / 0.2e1;
t217 = t258 - t263 / 0.2e1;
t236 = cos(pkin(6));
t197 = t213 * t209 + t217 * t232 + t236 * t268;
t198 = t213 * t208 + t236 * t210 - t217 * t235;
t202 = -t213 * t233 + t236 * t285;
t182 = t197 * t212 + t198 * t244 - t202 * t215;
t188 = -t197 * t283 + t202 * t284;
t169 = t182 * t238 - t188 * t243;
t154 = atan2(t288, t169);
t151 = sin(t154);
t152 = cos(t154);
t149 = t151 * t288 + t152 * t169;
t148 = 0.1e1 / t149 ^ 2;
t206 = t241 * t214 - t246 * t245;
t252 = t246 * t240 + t241 * t254;
t250 = -t206 * t232 + t209 * t252 - t241 * t259;
t251 = -t233 * t252 - t241 * t264;
t247 = t250 * t283 - t251 * t284;
t193 = -t206 * t235 - t208 * t252 + t241 * t273;
t248 = t193 * t244 - t212 * t250 + t215 * t251;
t163 = t238 * t248 - t243 * t247;
t282 = t148 * t163;
t281 = t152 * t288;
t164 = t238 * t247 + t243 * t248;
t175 = t193 * t239 + t211 * t251 + t216 * t250;
t237 = sin(qJ(6));
t242 = cos(qJ(6));
t158 = t164 * t242 + t175 * t237;
t156 = 0.1e1 / t158 ^ 2;
t157 = t164 * t237 - t175 * t242;
t280 = t156 * t157;
t168 = 0.1e1 / t169 ^ 2;
t279 = t288 * t168;
t278 = t163 ^ 2 * t148;
t277 = t175 * t243;
t274 = t233 * t215;
t267 = t157 ^ 2 * t156 + 0.1e1;
t265 = t233 * t284;
t255 = -t151 * t169 + t281;
t200 = t217 * t209 - t213 * t232;
t196 = t206 * t208 - t235 * t252;
t195 = t206 * t209 + t232 * t252;
t194 = t203 * t232 - t204 * t209;
t186 = -t195 * t283 - t206 * t265;
t181 = t197 * t216 - t198 * t239 + t202 * t211;
t180 = t195 * t212 + t196 * t244 + t206 * t274;
t179 = t206 * t233 * t211 - t195 * t216 + t196 * t239;
t178 = ((t217 * t208 + t213 * t235) * t244 + t200 * t212 + t217 * t274) * t238 - (-t200 * t283 - t217 * t265) * t243;
t170 = t182 * t243 + t188 * t238;
t167 = 0.1e1 / t169;
t166 = t180 * t243 + t186 * t238;
t165 = ((-t203 * t235 - t204 * t208) * t244 + t194 * t212 - t204 * t274) * t238 - (-t194 * t283 + t204 * t265) * t243;
t155 = 0.1e1 / t158;
t153 = 0.1e1 / (t168 * t288 ^ 2 + 0.1e1);
t150 = 0.1e1 / t267;
t147 = 0.1e1 / t149;
t146 = 0.1e1 / (0.1e1 + t278);
t145 = (-t167 * t287 - t181 * t279) * t238 * t153;
t144 = (-t165 * t167 - t178 * t279) * t153;
t143 = (t162 * t167 - t170 * t279) * t153;
t1 = [-t163 * t167 * t153, t144, 0, t145, t143, 0; (t288 * t147 - (-t151 + (-t167 * t281 + t151) * t153) * t278) * t146 ((t180 * t238 - t186 * t243) * t147 - (t144 * t255 - t151 * t165 + t152 * t178) * t282) * t146, 0 (-t175 * t238 * t147 - ((-t151 * t287 + t152 * t181) * t238 + t255 * t145) * t282) * t146 (t164 * t147 - (t143 * t255 + t151 * t162 + t152 * t170) * t282) * t146, 0; ((t162 * t237 - t242 * t287) * t155 - (t162 * t242 + t237 * t287) * t280) * t150 ((t166 * t237 - t179 * t242) * t155 - (t166 * t242 + t179 * t237) * t280) * t150, 0 ((-t237 * t277 - t242 * t248) * t155 - (t237 * t248 - t242 * t277) * t280) * t150 (-t237 * t155 + t242 * t280) * t163 * t150, t267 * t150;];
Ja_rot  = t1;
