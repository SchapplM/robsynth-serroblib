% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:33
% EndTime: 2019-02-26 22:06:34
% DurationCPUTime: 1.10s
% Computational Cost: add. (7177->123), mult. (11003->252), div. (691->12), fcn. (14028->11), ass. (0->106)
t210 = cos(pkin(6));
t212 = sin(qJ(1));
t213 = cos(qJ(2));
t269 = cos(qJ(1));
t234 = t269 * qJD(2);
t235 = t269 * qJD(1);
t209 = sin(pkin(6));
t240 = t209 * t269;
t211 = sin(qJ(2));
t253 = t212 * t211;
t241 = t210 * t253;
t251 = qJD(2) * t211;
t274 = -qJD(1) * t241 - t212 * t251 + (t210 * t234 + t235) * t213 - qJD(3) * t240;
t208 = qJ(3) + pkin(11);
t206 = sin(t208);
t207 = cos(t208);
t239 = t269 * t211;
t252 = t212 * t213;
t221 = -t210 * t239 - t252;
t177 = -t206 * t221 + t207 * t240;
t256 = t209 * t211;
t187 = t206 * t256 - t210 * t207;
t164 = atan2(-t177, t187);
t159 = sin(t164);
t160 = cos(t164);
t154 = -t159 * t177 + t160 * t187;
t152 = 0.1e1 / t154 ^ 2;
t238 = t269 * t213;
t195 = t238 - t241;
t255 = t209 * t212;
t225 = -t195 * t206 + t207 * t255;
t175 = t225 ^ 2;
t150 = t175 * t152 + 0.1e1;
t222 = -t210 * t252 - t239;
t169 = t221 * qJD(1) + t222 * qJD(2);
t183 = t195 * t207 + t206 * t255;
t229 = t209 * t235;
t155 = t183 * qJD(3) + t169 * t206 - t207 * t229;
t264 = t155 * t152;
t174 = t177 ^ 2;
t185 = 0.1e1 / t187 ^ 2;
t163 = t174 * t185 + 0.1e1;
t161 = 0.1e1 / t163;
t237 = qJD(1) * t255;
t250 = qJD(3) * t207;
t157 = t274 * t206 - t207 * t237 - t221 * t250;
t188 = t210 * t206 + t207 * t256;
t254 = t209 * t213;
t236 = qJD(2) * t254;
t172 = t188 * qJD(3) + t206 * t236;
t184 = 0.1e1 / t187;
t259 = t177 * t185;
t227 = -t157 * t184 + t172 * t259;
t143 = t227 * t161;
t228 = -t159 * t187 - t160 * t177;
t139 = t228 * t143 - t159 * t157 + t160 * t172;
t151 = 0.1e1 / t154;
t153 = t151 * t152;
t267 = t139 * t153;
t249 = 0.2e1 * (-t175 * t267 - t225 * t264) / t150 ^ 2;
t273 = t172 * t185;
t231 = t210 * t238;
t192 = t231 - t253;
t223 = -t184 * t192 + t254 * t259;
t272 = t206 * t223;
t158 = (qJD(3) * t221 + t237) * t206 + t274 * t207;
t189 = 0.1e1 / t222;
t190 = 0.1e1 / t222 ^ 2;
t271 = -0.2e1 * t177;
t270 = -0.2e1 * t225;
t261 = t184 * t273;
t266 = (t157 * t259 - t174 * t261) / t163 ^ 2;
t265 = t152 * t225;
t263 = t159 * t225;
t262 = t160 * t225;
t260 = t177 * t184;
t258 = t183 * t190;
t257 = t222 * t206;
t248 = -0.2e1 * t266;
t156 = t225 * qJD(3) + t169 * t207 + t206 * t229;
t176 = t183 ^ 2;
t167 = t176 * t190 + 0.1e1;
t168 = -qJD(1) * t231 - t213 * t234 + (qJD(2) * t210 + qJD(1)) * t253;
t191 = t189 * t190;
t247 = 0.2e1 * (-t176 * t191 * t168 + t156 * t258) / t167 ^ 2;
t246 = t153 * t270;
t245 = 0.2e1 * t183 * t191;
t244 = t184 * t266;
t243 = t152 * t263;
t242 = t152 * t262;
t233 = t261 * t271;
t179 = -t206 * t240 - t207 * t221;
t226 = -t179 * t184 + t188 * t259;
t220 = -t159 + (t160 * t260 + t159) * t161;
t173 = -t187 * qJD(3) + t207 * t236;
t170 = t222 * qJD(1) + t221 * qJD(2);
t165 = 0.1e1 / t167;
t148 = 0.1e1 / t150;
t147 = t161 * t272;
t144 = t226 * t161;
t142 = t220 * t225;
t141 = (-t159 * t192 + t160 * t254) * t206 + t228 * t147;
t140 = t228 * t144 - t159 * t179 + t160 * t188;
t138 = t226 * t248 + (t188 * t233 - t158 * t184 + (t157 * t188 + t172 * t179 + t173 * t177) * t185) * t161;
t136 = t248 * t272 + (t223 * t250 + (t233 * t254 - t170 * t184 + (t172 * t192 + (t157 * t213 - t177 * t251) * t209) * t185) * t206) * t161;
t1 = [t244 * t270 + (-t155 * t184 - t225 * t273) * t161, t136, t138, 0, 0, 0; t177 * t151 * t249 + (-t157 * t151 + (t139 * t177 + t142 * t155) * t152) * t148 - (-t142 * t152 * t249 + (-0.2e1 * t142 * t267 + (-t143 * t161 * t260 + t248) * t243 + (t244 * t271 - t143 + (t143 - t227) * t161) * t242 - t220 * t264) * t148) * t225 (-t141 * t265 - t151 * t257) * t249 + (-t141 * t264 + (t168 * t206 + t222 * t250) * t151 + (t141 * t246 - t152 * t257) * t139 + (-t136 * t177 - t147 * t157 + (-t206 * t251 + t213 * t250) * t209 + (-t147 * t187 - t192 * t206) * t143) * t242 + (-t192 * t250 - t136 * t187 - t147 * t172 - t170 * t206 + (t147 * t177 - t206 * t254) * t143) * t243) * t148 (-t140 * t265 - t151 * t183) * t249 + (t140 * t139 * t246 + t156 * t151 + (-t183 * t139 - t140 * t155 + (-t138 * t177 - t144 * t157 + t173 + (-t144 * t187 - t179) * t143) * t262 + (-t138 * t187 - t144 * t172 - t158 + (t144 * t177 - t188) * t143) * t263) * t152) * t148, 0, 0, 0; (-t179 * t189 + t192 * t258) * t247 + (t158 * t189 + t192 * t168 * t245 + (-t192 * t156 - t168 * t179 - t170 * t183) * t190) * t165 (t189 * t207 * t222 + t195 * t258) * t247 + (qJD(3) * t189 * t257 + (-t156 * t195 - t169 * t183) * t190 + (t195 * t245 + (t190 * t222 - t189) * t207) * t168) * t165, t225 * t189 * t247 + (t168 * t190 * t225 + t155 * t189) * t165, 0, 0, 0;];
JaD_rot  = t1;
