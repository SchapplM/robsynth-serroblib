% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:09
% EndTime: 2019-02-26 22:05:10
% DurationCPUTime: 1.06s
% Computational Cost: add. (2275->116), mult. (4276->247), div. (493->12), fcn. (5073->11), ass. (0->112)
t271 = qJD(3) - qJD(6);
t192 = sin(qJ(2));
t185 = t192 ^ 2;
t195 = cos(qJ(2));
t188 = 0.1e1 / t195 ^ 2;
t249 = t185 * t188;
t193 = sin(qJ(1));
t186 = t193 ^ 2;
t178 = t186 * t249 + 0.1e1;
t187 = 0.1e1 / t195;
t246 = t187 * t192;
t269 = t192 * t249;
t205 = qJD(2) * (t187 * t269 + t246);
t196 = cos(qJ(1));
t237 = qJD(1) * t196;
t247 = t185 * t193;
t213 = t237 * t247;
t252 = (t186 * t205 + t188 * t213) / t178 ^ 2;
t270 = -0.2e1 * t252;
t220 = 0.1e1 + t249;
t268 = t193 * t220;
t183 = qJ(3) + pkin(10);
t182 = cos(t183);
t239 = t196 * t182;
t181 = sin(t183);
t244 = t193 * t181;
t171 = t195 * t239 + t244;
t194 = cos(qJ(6));
t267 = t271 * t194;
t191 = sin(qJ(6));
t266 = t271 * t191;
t241 = t193 * t195;
t206 = t181 * t241 + t239;
t234 = qJD(2) * t196;
t221 = t192 * t234;
t141 = t206 * qJD(1) - t171 * qJD(3) + t181 * t221;
t215 = -qJD(1) * t195 + qJD(3);
t216 = qJD(3) * t195 - qJD(1);
t240 = t196 * t181;
t142 = -t216 * t240 + (t215 * t193 - t221) * t182;
t243 = t193 * t182;
t170 = t195 * t240 - t243;
t209 = t170 * t194 - t171 * t191;
t135 = t209 * qJD(6) - t141 * t191 + t142 * t194;
t155 = t170 * t191 + t171 * t194;
t147 = 0.1e1 / t155;
t207 = t181 * t191 + t182 * t194;
t208 = t181 * t194 - t182 * t191;
t148 = 0.1e1 / t155 ^ 2;
t256 = t148 * t209;
t265 = t208 * t147 - t207 * t256;
t242 = t193 * t192;
t177 = atan2(t242, t195);
t174 = cos(t177);
t173 = sin(t177);
t225 = t173 * t242;
t162 = t174 * t195 + t225;
t159 = 0.1e1 / t162;
t160 = 0.1e1 / t162 ^ 2;
t264 = 0.2e1 * t192;
t175 = 0.1e1 / t178;
t263 = t175 - 0.1e1;
t134 = t155 * qJD(6) + t141 * t194 + t142 * t191;
t146 = t209 ^ 2;
t139 = t146 * t148 + 0.1e1;
t149 = t147 * t148;
t259 = t135 * t149;
t262 = (-t134 * t256 - t146 * t259) / t139 ^ 2;
t190 = t196 ^ 2;
t248 = t185 * t190;
t158 = t160 * t248 + 0.1e1;
t235 = qJD(2) * t195;
t222 = t192 * t237;
t236 = qJD(2) * t193;
t145 = ((t193 * t235 + t222) * t187 + t236 * t249) * t175;
t250 = t174 * t192;
t132 = (t145 * t193 - qJD(2)) * t250 + (t222 + (-t145 + t236) * t195) * t173;
t260 = t132 * t159 * t160;
t261 = (-t248 * t260 + (t190 * t192 * t235 - t213) * t160) / t158 ^ 2;
t258 = t145 * t173;
t257 = t145 * t192;
t245 = t192 * t196;
t166 = t207 * t245;
t255 = t148 * t166;
t254 = t160 * t192;
t253 = t160 * t196;
t164 = t175 * t268;
t251 = t164 * t193;
t238 = qJD(1) * t193;
t229 = 0.2e1 * t262;
t228 = -0.2e1 * t260;
t227 = -0.2e1 * t149 * t209;
t226 = t160 * t245;
t224 = t175 * t185 * t187;
t219 = -0.2e1 * t192 * t261;
t218 = t135 * t227;
t217 = t187 * t270;
t214 = t193 * t224;
t212 = t220 * t196;
t169 = -t182 * t241 + t240;
t210 = -t169 * t191 - t194 * t206;
t151 = t169 * t194 - t191 * t206;
t204 = t192 * t236 + t215 * t196;
t165 = t208 * t245;
t156 = 0.1e1 / t158;
t144 = t204 * t182 + t216 * t244;
t143 = t204 * t181 - t216 * t243;
t140 = (-t263 * t192 * t173 + t174 * t214) * t196;
t137 = 0.1e1 / t139;
t136 = t173 * t241 - t250 + (-t173 * t195 + t174 * t242) * t164;
t133 = t268 * t270 + (qJD(1) * t212 + 0.2e1 * t193 * t205) * t175;
t1 = [t217 * t245 + (qJD(2) * t212 - t238 * t246) * t175, t133, 0, 0, 0, 0; (t159 * t219 + (t159 * t235 + (-qJD(1) * t140 - t132) * t254) * t156) * t193 + (t160 * t219 * t140 + (((-t145 * t214 - t263 * t235 + t252 * t264) * t173 + (t217 * t247 + t257 + (-t257 + (t264 + t269) * t236) * t175) * t174) * t226 + (t160 * t235 + t192 * t228) * t140 + (t159 + ((-t186 + t190) * t174 * t224 + t263 * t225) * t160) * t192 * qJD(1)) * t156) * t196, 0.2e1 * (-t136 * t254 + t159 * t195) * t196 * t261 + ((t159 * t238 + (qJD(2) * t136 + t132) * t253) * t195 + (t159 * t234 + (t133 * t174 * t193 - t173 * t236 - t251 * t258 + t258 + (qJD(2) * t173 + t174 * t237) * t164) * t226 + (-t160 * t238 + t196 * t228) * t136 + ((-t133 + t237) * t173 + ((-0.1e1 + t251) * qJD(2) + (-t164 + t193) * t145) * t174) * t195 * t253) * t192) * t156, 0, 0, 0, 0; (t147 * t210 - t151 * t256) * t229 + ((t151 * qJD(6) - t143 * t194 + t144 * t191) * t147 + t151 * t218 + (t210 * t135 + (t210 * qJD(6) + t143 * t191 + t144 * t194) * t209 - t151 * t134) * t148) * t137 (-t147 * t165 + t209 * t255) * t229 + (t134 * t255 + (-t165 * t148 - t166 * t227) * t135 + t265 * t195 * t234 + (-t265 * t238 + ((t267 * t147 - t266 * t256) * t182 + (t266 * t147 + t267 * t256) * t181) * t196) * t192) * t137 (t147 * t155 + t209 * t256) * t229 + (-t135 * t147 - t209 * t218 + (0.2e1 * t209 * t134 + t155 * t135) * t148) * t137, 0, 0, -0.2e1 * t262 - 0.2e1 * (t134 * t148 * t137 - (-t137 * t259 - t148 * t262) * t209) * t209;];
JaD_rot  = t1;
