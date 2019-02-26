% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR10
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
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:23
% EndTime: 2019-02-26 22:21:24
% DurationCPUTime: 0.95s
% Computational Cost: add. (2039->120), mult. (4806->252), div. (511->12), fcn. (5689->11), ass. (0->114)
t208 = sin(qJ(2));
t200 = t208 ^ 2;
t211 = cos(qJ(2));
t203 = 0.1e1 / t211 ^ 2;
t264 = t200 * t203;
t209 = sin(qJ(1));
t201 = t209 ^ 2;
t192 = t201 * t264 + 0.1e1;
t202 = 0.1e1 / t211;
t261 = t202 * t208;
t287 = t208 * t264;
t220 = qJD(2) * (t202 * t287 + t261);
t212 = cos(qJ(1));
t253 = qJD(1) * t212;
t262 = t200 * t209;
t227 = t253 * t262;
t269 = (t201 * t220 + t203 * t227) / t192 ^ 2;
t288 = -0.2e1 * t269;
t238 = 0.1e1 + t264;
t286 = t209 * t238;
t198 = qJD(5) + qJD(6);
t210 = cos(qJ(3));
t285 = (qJD(3) - t198) * t210;
t207 = sin(qJ(3));
t249 = qJD(3) * t207;
t284 = -t198 * t207 + t249;
t256 = t212 * t207;
t258 = t209 * t210;
t185 = t211 * t256 - t258;
t255 = t212 * t210;
t241 = t211 * t255;
t186 = t209 * t207 + t241;
t206 = qJ(5) + qJ(6);
t196 = sin(t206);
t197 = cos(t206);
t170 = t185 * t196 + t186 * t197;
t162 = 0.1e1 / t170;
t223 = t196 * t207 + t197 * t210;
t224 = t196 * t210 - t197 * t207;
t163 = 0.1e1 / t170 ^ 2;
t225 = -t185 * t197 + t186 * t196;
t273 = t163 * t225;
t283 = t224 * t162 - t223 * t273;
t229 = -qJD(1) * t211 + qJD(3);
t230 = qJD(3) * t211 - qJD(1);
t250 = qJD(2) * t212;
t239 = t208 * t250;
t233 = t185 * t198 - t230 * t256 + (t229 * t209 - t239) * t210;
t257 = t209 * t211;
t221 = t207 * t257 + t255;
t234 = t221 * qJD(1) - qJD(3) * t241 + t186 * t198 + t207 * t239 - t209 * t249;
t282 = t234 * t196 - t233 * t197;
t259 = t209 * t208;
t191 = atan2(t259, t211);
t188 = cos(t191);
t187 = sin(t191);
t243 = t187 * t259;
t177 = t188 * t211 + t243;
t174 = 0.1e1 / t177;
t175 = 0.1e1 / t177 ^ 2;
t281 = 0.2e1 * t208;
t189 = 0.1e1 / t192;
t280 = t189 - 0.1e1;
t149 = t233 * t196 + t234 * t197;
t161 = t225 ^ 2;
t154 = t161 * t163 + 0.1e1;
t164 = t162 * t163;
t276 = t282 * t164;
t279 = (t149 * t273 + t161 * t276) / t154 ^ 2;
t205 = t212 ^ 2;
t263 = t200 * t205;
t173 = t175 * t263 + 0.1e1;
t251 = qJD(2) * t211;
t240 = t208 * t253;
t252 = qJD(2) * t209;
t156 = ((t209 * t251 + t240) * t202 + t252 * t264) * t189;
t267 = t188 * t208;
t147 = (t156 * t209 - qJD(2)) * t267 + (t240 + (-t156 + t252) * t211) * t187;
t277 = t147 * t174 * t175;
t278 = (-t263 * t277 + (t205 * t208 * t251 - t227) * t175) / t173 ^ 2;
t275 = t156 * t187;
t274 = t156 * t208;
t260 = t208 * t212;
t181 = t223 * t260;
t272 = t163 * t181;
t271 = t175 * t208;
t270 = t175 * t212;
t179 = t189 * t286;
t268 = t179 * t209;
t254 = qJD(1) * t209;
t247 = 0.2e1 * t279;
t246 = -0.2e1 * t277;
t245 = 0.2e1 * t164 * t225;
t244 = t175 * t260;
t242 = t189 * t200 * t202;
t237 = -0.2e1 * t208 * t278;
t236 = t282 * t245;
t235 = t202 * t288;
t184 = -t210 * t257 + t256;
t222 = t229 * t212;
t232 = t184 * t198 + t230 * t258 - (t208 * t252 + t222) * t207;
t231 = -t221 * t198 + t210 * t222 + (qJD(2) * t208 * t210 + t230 * t207) * t209;
t228 = t209 * t242;
t226 = t238 * t212;
t180 = t224 * t260;
t171 = 0.1e1 / t173;
t166 = t184 * t197 - t196 * t221;
t165 = t184 * t196 + t197 * t221;
t155 = (-t280 * t208 * t187 + t188 * t228) * t212;
t152 = 0.1e1 / t154;
t151 = t187 * t257 - t267 + (-t187 * t211 + t188 * t259) * t179;
t148 = t286 * t288 + (qJD(1) * t226 + 0.2e1 * t209 * t220) * t189;
t144 = -0.2e1 * t279 + 0.2e1 * (t149 * t163 * t152 + (t152 * t276 - t163 * t279) * t225) * t225;
t1 = [t235 * t260 + (qJD(2) * t226 - t254 * t261) * t189, t148, 0, 0, 0, 0; (t174 * t237 + (t174 * t251 + (-qJD(1) * t155 - t147) * t271) * t171) * t209 + (t175 * t237 * t155 + (((-t156 * t228 - t280 * t251 + t269 * t281) * t187 + (t235 * t262 + t274 + (-t274 + (t281 + t287) * t252) * t189) * t188) * t244 + (t175 * t251 + t208 * t246) * t155 + (t174 + ((-t201 + t205) * t188 * t242 + t280 * t243) * t175) * t208 * qJD(1)) * t171) * t212, 0.2e1 * (-t151 * t271 + t174 * t211) * t212 * t278 + ((t174 * t254 + (qJD(2) * t151 + t147) * t270) * t211 + (t174 * t250 + (t148 * t188 * t209 - t187 * t252 - t268 * t275 + t275 + (qJD(2) * t187 + t188 * t253) * t179) * t244 + (-t175 * t254 + t212 * t246) * t151 + ((-t148 + t253) * t187 + ((-0.1e1 + t268) * qJD(2) + (-t179 + t209) * t156) * t188) * t211 * t270) * t208) * t171, 0, 0, 0, 0; (-t162 * t165 + t166 * t273) * t247 + ((t231 * t196 + t232 * t197) * t162 - t166 * t236 + (t165 * t282 - (-t232 * t196 + t231 * t197) * t225 - t166 * t149) * t163) * t152 (t162 * t180 - t225 * t272) * t247 + (t149 * t272 - (t180 * t163 - t181 * t245) * t282 - t283 * t211 * t250 + (t283 * t254 + ((t285 * t162 - t284 * t273) * t197 + (t284 * t162 + t285 * t273) * t196) * t212) * t208) * t152 (t162 * t170 + t225 * t273) * t247 + (t282 * t162 - t225 * t236 + (-0.2e1 * t225 * t149 - t170 * t282) * t163) * t152, 0, t144, t144;];
JaD_rot  = t1;
