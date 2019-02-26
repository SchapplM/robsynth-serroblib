% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:36
% EndTime: 2019-02-26 22:10:37
% DurationCPUTime: 1.17s
% Computational Cost: add. (9575->126), mult. (8382->275), div. (1515->15), fcn. (10508->9), ass. (0->118)
t198 = pkin(10) + qJ(5);
t194 = sin(t198);
t201 = qJ(2) + qJ(3);
t197 = cos(t201);
t195 = cos(t198);
t202 = cos(qJ(1));
t251 = t202 * t195;
t273 = sin(qJ(1));
t176 = t273 * t194 + t197 * t251;
t170 = 0.1e1 / t176 ^ 2;
t196 = sin(t201);
t190 = t196 ^ 2;
t200 = t202 ^ 2;
t258 = t190 * t200;
t238 = t170 * t258;
t166 = 0.1e1 + t238;
t226 = qJD(1) * t273;
t199 = qJD(2) + qJD(3);
t254 = t199 * t202;
t232 = t196 * t254;
t212 = t197 * t226 + t232;
t225 = t273 * qJD(5);
t252 = t202 * t194;
t155 = (-qJD(5) * t197 + qJD(1)) * t252 + (t225 - t212) * t195;
t169 = 0.1e1 / t176;
t268 = t155 * t169 * t170;
t220 = t258 * t268;
t233 = t196 * t199 * t200;
t276 = (-t220 + (-t190 * t202 * t226 + t197 * t233) * t170) / t166 ^ 2;
t256 = t196 * t202;
t230 = t273 * t197;
t172 = t194 * t230 + t251;
t217 = t194 * t225;
t247 = qJD(5) * t202;
t228 = t195 * t247;
t154 = t172 * qJD(1) + t194 * t232 - t197 * t228 - t217;
t175 = -t273 * t195 + t197 * t252;
t187 = 0.1e1 / t194;
t188 = 0.1e1 / t194 ^ 2;
t191 = 0.1e1 / t196;
t192 = 0.1e1 / t196 ^ 2;
t255 = t197 * t199;
t234 = t192 * t255;
t249 = qJD(5) * t195;
t261 = t187 * t191;
t275 = (t188 * t191 * t249 + t187 * t234) * t175 + t154 * t261;
t257 = t196 * t194;
t162 = atan2(-t172, t257);
t159 = cos(t162);
t158 = sin(t162);
t267 = t158 * t172;
t153 = t159 * t257 - t267;
t150 = 0.1e1 / t153;
t151 = 0.1e1 / t153 ^ 2;
t274 = 0.2e1 * t175;
t167 = t172 ^ 2;
t260 = t188 * t192;
t163 = t167 * t260 + 0.1e1;
t160 = 0.1e1 / t163;
t248 = qJD(5) * t196;
t213 = t194 * t255 + t195 * t248;
t236 = t172 * t260;
t218 = t195 * t226;
t231 = t273 * t196;
t219 = t199 * t231;
t250 = qJD(1) * t202;
t156 = t195 * t225 * t197 - t218 + (t250 * t197 - t219 - t247) * t194;
t239 = t156 * t261;
t142 = (t213 * t236 - t239) * t160;
t210 = -t142 * t172 + t213;
t138 = (-t142 * t257 - t156) * t158 + t210 * t159;
t152 = t150 * t151;
t272 = t138 * t152;
t189 = t187 * t188;
t193 = t191 / t190;
t229 = t192 * t249;
t271 = (t156 * t236 + (-t188 * t193 * t255 - t189 * t229) * t167) / t163 ^ 2;
t270 = t151 * t175;
t269 = t154 * t151;
t266 = t158 * t175;
t265 = t158 * t196;
t264 = t159 * t172;
t263 = t159 * t175;
t262 = t159 * t197;
t259 = t188 * t195;
t253 = t202 * t150;
t168 = t175 ^ 2;
t148 = t151 * t168 + 0.1e1;
t246 = 0.2e1 * (-t168 * t272 - t175 * t269) / t148 ^ 2;
t245 = -0.2e1 * t271;
t244 = 0.2e1 * t276;
t243 = t152 * t274;
t242 = t191 * t271;
t241 = t151 * t266;
t237 = t172 * t261;
t235 = t187 * t192 * t197;
t215 = t172 * t235 + t273;
t149 = t215 * t160;
t227 = t273 - t149;
t224 = t150 * t246;
t223 = t151 * t246;
t222 = t256 * t274;
t221 = t187 * t242;
t174 = t195 * t230 - t252;
t216 = t172 * t259 - t174 * t187;
t214 = t170 * t174 * t202 - t273 * t169;
t164 = 0.1e1 / t166;
t157 = t176 * qJD(1) - t195 * t219 - t197 * t217 - t228;
t146 = 0.1e1 / t148;
t145 = t216 * t191 * t160;
t141 = (-t158 + (t159 * t237 + t158) * t160) * t175;
t140 = -t149 * t264 + (t227 * t265 + t262) * t194;
t139 = t159 * t195 * t196 - t158 * t174 + (-t158 * t257 - t264) * t145;
t137 = t215 * t245 + (t156 * t235 + t250 + (-t188 * t197 * t229 + (-0.2e1 * t193 * t197 ^ 2 - t191) * t199 * t187) * t172) * t160;
t135 = (t169 * t197 * t202 + t195 * t238) * t244 + (0.2e1 * t195 * t220 + t212 * t169 + ((t155 * t202 - 0.2e1 * t195 * t233) * t197 + (qJD(5) * t194 * t200 + 0.2e1 * t202 * t218) * t190) * t170) * t164;
t134 = -0.2e1 * t216 * t242 + (-t216 * t234 + (t156 * t259 - t157 * t187 + (t174 * t259 + (-0.2e1 * t189 * t195 ^ 2 - t187) * t172) * qJD(5)) * t191) * t160;
t133 = t140 * t175 * t223 + (-(-t137 * t264 + (t142 * t267 - t156 * t159) * t149) * t270 + (t138 * t243 + t269) * t140 + (-t196 * t253 - (-t149 * t265 + t158 * t231 + t262) * t270) * t249) * t146 + (t224 * t256 + ((-t199 * t253 - (t227 * t199 - t142) * t241) * t197 + (t150 * t226 + (t202 * t138 - (-t137 + t250) * t266 - (t227 * t142 - t199) * t263) * t151) * t196) * t146) * t194;
t1 = [t275 * t160 + t221 * t274, t137, t137, 0, t134, 0; t172 * t224 + (-t156 * t150 + (t138 * t172 + t141 * t154) * t151) * t146 + (t141 * t223 + (0.2e1 * t141 * t272 + (t154 * t160 - t154 - (-t142 * t160 * t237 + t245) * t175) * t151 * t158 + (-(-0.2e1 * t172 * t221 - t142) * t270 + (-(t142 + t239) * t175 + t275 * t172) * t151 * t160) * t159) * t146) * t175, t133, t133, 0 (t139 * t270 - t150 * t176) * t246 + (t139 * t269 + t155 * t150 + (t139 * t243 - t151 * t176) * t138 - (-t194 * t248 + t195 * t255 - t134 * t172 - t145 * t156 + (-t145 * t257 - t174) * t142) * t151 * t263 - (-t157 + (-t134 * t194 - t142 * t195) * t196 - t210 * t145) * t241) * t146, 0; t214 * t196 * t244 + (-t214 * t255 + ((qJD(1) * t169 + 0.2e1 * t174 * t268) * t202 + (-t273 * t155 - t157 * t202 + t174 * t226) * t170) * t196) * t164, t135, t135, 0, t170 * t222 * t276 + (t222 * t268 + (t154 * t256 + (t196 * t226 - t197 * t254) * t175) * t170) * t164, 0;];
JaD_rot  = t1;
