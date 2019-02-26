% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:57
% EndTime: 2019-02-26 22:52:59
% DurationCPUTime: 1.74s
% Computational Cost: add. (6679->105), mult. (19296->240), div. (145->9), fcn. (26142->21), ass. (0->114)
t227 = cos(pkin(6));
t232 = sin(qJ(2));
t284 = sin(qJ(1));
t256 = t284 * t232;
t237 = cos(qJ(2));
t238 = cos(qJ(1));
t259 = t238 * t237;
t216 = -t227 * t259 + t256;
t255 = t284 * t237;
t260 = t238 * t232;
t217 = t227 * t260 + t255;
t231 = sin(qJ(3));
t236 = cos(qJ(3));
t223 = sin(pkin(7));
t224 = sin(pkin(6));
t269 = t224 * t238;
t258 = t223 * t269;
t226 = cos(pkin(7));
t267 = t226 * t231;
t205 = t216 * t267 - t217 * t236 + t231 * t258;
t230 = sin(qJ(4));
t235 = cos(qJ(4));
t204 = (t216 * t226 + t258) * t236 + t217 * t231;
t213 = -t216 * t223 + t226 * t269;
t222 = sin(pkin(8));
t225 = cos(pkin(8));
t252 = t204 * t225 + t213 * t222;
t186 = t205 * t235 + t252 * t230;
t229 = sin(qJ(5));
t234 = cos(qJ(5));
t241 = t204 * t222 - t213 * t225;
t290 = t186 * t229 + t241 * t234;
t172 = t186 * t234 - t241 * t229;
t287 = t205 * t230 - t252 * t235;
t262 = t232 * t236;
t264 = t231 * t237;
t270 = t223 * t227;
t212 = t231 * t270 + (t226 * t264 + t262) * t224;
t261 = t236 * t237;
t265 = t231 * t232;
t211 = t236 * t270 + (t226 * t261 - t265) * t224;
t215 = -t224 * t237 * t223 + t227 * t226;
t251 = t211 * t225 + t215 * t222;
t196 = t212 * t235 + t251 * t230;
t202 = -t211 * t222 + t215 * t225;
t181 = t196 * t229 - t202 * t234;
t168 = atan2(t290, t181);
t161 = sin(t168);
t162 = cos(t168);
t159 = t161 * t290 + t162 * t181;
t158 = 0.1e1 / t159 ^ 2;
t218 = -t227 * t256 + t259;
t248 = t227 * t255 + t260;
t257 = t224 * t284;
t245 = t223 * t257 - t248 * t226;
t243 = t218 * t231 - t245 * t236;
t246 = t248 * t223 + t226 * t257;
t244 = t246 * t222;
t206 = t218 * t236 + t245 * t231;
t277 = t206 * t235;
t188 = t277 + (-t243 * t225 + t244) * t230;
t239 = t243 * t222 + t246 * t225;
t173 = t188 * t229 - t239 * t234;
t283 = t158 * t173;
t282 = t162 * t290;
t174 = t188 * t234 + t239 * t229;
t242 = t243 * t235;
t187 = t206 * t230 + t225 * t242 - t235 * t244;
t228 = sin(qJ(6));
t233 = cos(qJ(6));
t167 = t174 * t233 + t187 * t228;
t164 = 0.1e1 / t167 ^ 2;
t166 = t174 * t228 - t187 * t233;
t281 = t164 * t166;
t180 = 0.1e1 / t181 ^ 2;
t280 = t290 * t180;
t279 = t173 ^ 2 * t158;
t278 = t187 * t234;
t273 = t222 * t234;
t272 = t223 * t222;
t271 = t223 * t225;
t268 = t225 * t230;
t266 = t226 * t236;
t263 = t232 * t223;
t254 = t166 ^ 2 * t164 + 0.1e1;
t253 = -t161 * t181 + t282;
t208 = -t218 * t266 + t248 * t231;
t250 = t208 * t225 + t218 * t272;
t209 = -t218 * t267 - t248 * t236;
t207 = t216 * t231 - t217 * t266;
t200 = -t208 * t222 + t218 * t271;
t195 = -t212 * t230 + t251 * t235;
t194 = (t211 * t235 - t212 * t268) * t229 - t212 * t273;
t193 = ((t229 * t268 + t273) * (-t226 * t262 - t264) + ((-t226 * t265 + t261) * t235 + t222 * t230 * t263) * t229 - t225 * t234 * t263) * t224;
t192 = -t206 * t268 - t242;
t191 = t225 * t277 - t243 * t230;
t190 = t209 * t235 + t250 * t230;
t189 = t209 * t230 - t250 * t235;
t182 = t196 * t234 + t202 * t229;
t179 = 0.1e1 / t181;
t178 = t206 * t222 * t229 + t192 * t234;
t177 = (-t204 * t235 + t205 * t268) * t229 + t205 * t273;
t176 = t190 * t234 + t200 * t229;
t175 = ((-t216 * t236 - t217 * t267) * t235 + (t207 * t225 + t217 * t272) * t230) * t229 - (-t207 * t222 + t217 * t271) * t234;
t165 = 0.1e1 / (t180 * t290 ^ 2 + 0.1e1);
t163 = 0.1e1 / t167;
t160 = 0.1e1 / t254;
t157 = 0.1e1 / t159;
t156 = 0.1e1 / (0.1e1 + t279);
t155 = (-t179 * t287 - t195 * t280) * t229 * t165;
t154 = (-t177 * t179 - t194 * t280) * t165;
t153 = (-t175 * t179 - t193 * t280) * t165;
t152 = (t172 * t179 - t182 * t280) * t165;
t1 = [-t173 * t179 * t165, t153, t154, t155, t152, 0; (t290 * t157 - (-t161 + (-t179 * t282 + t161) * t165) * t279) * t156 ((t190 * t229 - t200 * t234) * t157 - (t253 * t153 - t161 * t175 + t162 * t193) * t283) * t156 ((t192 * t229 - t206 * t273) * t157 - (t253 * t154 - t161 * t177 + t162 * t194) * t283) * t156 (-t187 * t229 * t157 - ((-t161 * t287 + t162 * t195) * t229 + t253 * t155) * t283) * t156 (t174 * t157 - (t253 * t152 + t161 * t172 + t162 * t182) * t283) * t156, 0; ((t172 * t228 - t233 * t287) * t163 - (t172 * t233 + t228 * t287) * t281) * t160 ((t176 * t228 - t189 * t233) * t163 - (t176 * t233 + t189 * t228) * t281) * t160 ((t178 * t228 - t191 * t233) * t163 - (t178 * t233 + t191 * t228) * t281) * t160 ((-t188 * t233 - t228 * t278) * t163 - (t188 * t228 - t233 * t278) * t281) * t160 (-t228 * t163 + t233 * t281) * t173 * t160, t254 * t160;];
Ja_rot  = t1;
