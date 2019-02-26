% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:53
% EndTime: 2019-02-26 22:24:54
% DurationCPUTime: 1.09s
% Computational Cost: add. (9575->126), mult. (8382->275), div. (1515->15), fcn. (10508->9), ass. (0->118)
t202 = qJ(4) + pkin(10);
t197 = sin(t202);
t204 = qJ(2) + qJ(3);
t200 = cos(t204);
t198 = cos(t202);
t205 = cos(qJ(1));
t254 = t205 * t198;
t276 = sin(qJ(1));
t179 = t276 * t197 + t200 * t254;
t173 = 0.1e1 / t179 ^ 2;
t199 = sin(t204);
t193 = t199 ^ 2;
t203 = t205 ^ 2;
t261 = t193 * t203;
t241 = t173 * t261;
t169 = 0.1e1 + t241;
t229 = qJD(1) * t276;
t201 = qJD(2) + qJD(3);
t257 = t201 * t205;
t235 = t199 * t257;
t215 = t200 * t229 + t235;
t228 = t276 * qJD(4);
t255 = t205 * t197;
t158 = (-qJD(4) * t200 + qJD(1)) * t255 + (t228 - t215) * t198;
t172 = 0.1e1 / t179;
t271 = t158 * t172 * t173;
t223 = t261 * t271;
t236 = t199 * t201 * t203;
t279 = (-t223 + (-t193 * t205 * t229 + t200 * t236) * t173) / t169 ^ 2;
t259 = t199 * t205;
t233 = t276 * t200;
t175 = t197 * t233 + t254;
t220 = t197 * t228;
t250 = qJD(4) * t205;
t231 = t198 * t250;
t157 = t175 * qJD(1) + t197 * t235 - t200 * t231 - t220;
t178 = -t276 * t198 + t200 * t255;
t190 = 0.1e1 / t197;
t191 = 0.1e1 / t197 ^ 2;
t194 = 0.1e1 / t199;
t195 = 0.1e1 / t199 ^ 2;
t258 = t200 * t201;
t237 = t195 * t258;
t252 = qJD(4) * t198;
t264 = t190 * t194;
t278 = (t191 * t194 * t252 + t190 * t237) * t178 + t157 * t264;
t260 = t199 * t197;
t165 = atan2(-t175, t260);
t162 = cos(t165);
t161 = sin(t165);
t270 = t161 * t175;
t156 = t162 * t260 - t270;
t153 = 0.1e1 / t156;
t154 = 0.1e1 / t156 ^ 2;
t277 = 0.2e1 * t178;
t170 = t175 ^ 2;
t263 = t191 * t195;
t166 = t170 * t263 + 0.1e1;
t163 = 0.1e1 / t166;
t251 = qJD(4) * t199;
t216 = t197 * t258 + t198 * t251;
t239 = t175 * t263;
t221 = t198 * t229;
t234 = t276 * t199;
t222 = t201 * t234;
t253 = qJD(1) * t205;
t159 = t198 * t228 * t200 - t221 + (t253 * t200 - t222 - t250) * t197;
t242 = t159 * t264;
t145 = (t216 * t239 - t242) * t163;
t213 = -t145 * t175 + t216;
t141 = (-t145 * t260 - t159) * t161 + t213 * t162;
t155 = t153 * t154;
t275 = t141 * t155;
t192 = t190 * t191;
t196 = t194 / t193;
t232 = t195 * t252;
t274 = (t159 * t239 + (-t191 * t196 * t258 - t192 * t232) * t170) / t166 ^ 2;
t273 = t154 * t178;
t272 = t157 * t154;
t269 = t161 * t178;
t268 = t161 * t199;
t267 = t162 * t175;
t266 = t162 * t178;
t265 = t162 * t200;
t262 = t191 * t198;
t256 = t205 * t153;
t171 = t178 ^ 2;
t151 = t154 * t171 + 0.1e1;
t249 = 0.2e1 * (-t171 * t275 - t178 * t272) / t151 ^ 2;
t248 = -0.2e1 * t274;
t247 = 0.2e1 * t279;
t246 = t155 * t277;
t245 = t194 * t274;
t244 = t154 * t269;
t240 = t175 * t264;
t238 = t190 * t195 * t200;
t218 = t175 * t238 + t276;
t152 = t218 * t163;
t230 = t276 - t152;
t227 = t153 * t249;
t226 = t154 * t249;
t225 = t259 * t277;
t224 = t190 * t245;
t177 = t198 * t233 - t255;
t219 = t175 * t262 - t177 * t190;
t217 = t173 * t177 * t205 - t276 * t172;
t167 = 0.1e1 / t169;
t160 = t179 * qJD(1) - t198 * t222 - t200 * t220 - t231;
t149 = 0.1e1 / t151;
t148 = t219 * t194 * t163;
t144 = (-t161 + (t162 * t240 + t161) * t163) * t178;
t143 = -t152 * t267 + (t230 * t268 + t265) * t197;
t142 = t162 * t198 * t199 - t161 * t177 + (-t161 * t260 - t267) * t148;
t140 = t218 * t248 + (t159 * t238 + t253 + (-t191 * t200 * t232 + (-0.2e1 * t196 * t200 ^ 2 - t194) * t201 * t190) * t175) * t163;
t138 = (t172 * t200 * t205 + t198 * t241) * t247 + (0.2e1 * t198 * t223 + t215 * t172 + ((t158 * t205 - 0.2e1 * t198 * t236) * t200 + (qJD(4) * t197 * t203 + 0.2e1 * t205 * t221) * t193) * t173) * t167;
t137 = -0.2e1 * t219 * t245 + (-t219 * t237 + (t159 * t262 - t160 * t190 + (t177 * t262 + (-0.2e1 * t192 * t198 ^ 2 - t190) * t175) * qJD(4)) * t194) * t163;
t136 = t143 * t178 * t226 + (-(-t140 * t267 + (t145 * t270 - t159 * t162) * t152) * t273 + (t141 * t246 + t272) * t143 + (-t199 * t256 - (-t152 * t268 + t161 * t234 + t265) * t273) * t252) * t149 + (t227 * t259 + ((-t201 * t256 - (t230 * t201 - t145) * t244) * t200 + (t153 * t229 + (t205 * t141 - (-t140 + t253) * t269 - (t230 * t145 - t201) * t266) * t154) * t199) * t149) * t197;
t1 = [t163 * t278 + t224 * t277, t140, t140, t137, 0, 0; t175 * t227 + (-t159 * t153 + (t141 * t175 + t144 * t157) * t154) * t149 + (t144 * t226 + (0.2e1 * t144 * t275 + (t157 * t163 - t157 - (-t145 * t163 * t240 + t248) * t178) * t154 * t161 + (-(-0.2e1 * t175 * t224 - t145) * t273 + (-(t145 + t242) * t178 + t278 * t175) * t154 * t163) * t162) * t149) * t178, t136, t136 (t142 * t273 - t153 * t179) * t249 + (t142 * t272 + t158 * t153 + (t142 * t246 - t179 * t154) * t141 - (-t197 * t251 + t198 * t258 - t137 * t175 - t148 * t159 + (-t148 * t260 - t177) * t145) * t154 * t266 - (-t160 + (-t137 * t197 - t145 * t198) * t199 - t213 * t148) * t244) * t149, 0, 0; t217 * t199 * t247 + (-t217 * t258 + ((qJD(1) * t172 + 0.2e1 * t177 * t271) * t205 + (-t276 * t158 - t160 * t205 + t177 * t229) * t173) * t199) * t167, t138, t138, t173 * t225 * t279 + (t225 * t271 + (t157 * t259 + (t199 * t229 - t200 * t257) * t178) * t173) * t167, 0, 0;];
JaD_rot  = t1;
