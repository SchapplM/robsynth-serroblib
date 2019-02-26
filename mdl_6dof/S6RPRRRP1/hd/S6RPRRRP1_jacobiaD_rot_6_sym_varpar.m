% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:51
% EndTime: 2019-02-26 21:07:52
% DurationCPUTime: 1.11s
% Computational Cost: add. (8522->127), mult. (8382->276), div. (1515->15), fcn. (10508->9), ass. (0->119)
t198 = qJ(3) + qJ(4);
t193 = cos(t198);
t199 = sin(qJ(5));
t247 = qJ(1) + pkin(10);
t229 = sin(t247);
t222 = t229 * t199;
t191 = cos(t247);
t200 = cos(qJ(5));
t257 = t191 * t200;
t175 = t193 * t257 + t222;
t169 = 0.1e1 / t175 ^ 2;
t186 = t191 ^ 2;
t192 = sin(t198);
t187 = t192 ^ 2;
t262 = t186 * t187;
t236 = t169 * t262;
t159 = 0.1e1 + t236;
t219 = qJD(1) * t229;
t215 = t193 * t219;
t218 = t229 * qJD(5);
t194 = qJD(3) + qJD(4);
t252 = t194 * t200;
t154 = (-t215 + t218) * t200 + (-t192 * t252 + (-qJD(5) * t193 + qJD(1)) * t199) * t191;
t168 = 0.1e1 / t175;
t269 = t154 * t168 * t169;
t224 = t262 * t269;
t255 = t193 * t194;
t276 = (-t224 + (t186 * t192 * t255 - t187 * t191 * t219) * t169) / t159 ^ 2;
t259 = t191 * t192;
t171 = t193 * t222 + t257;
t214 = t199 * t218;
t248 = qJD(5) * t200;
t230 = t191 * t248;
t258 = t191 * t194;
t234 = t192 * t258;
t153 = qJD(1) * t171 - t193 * t230 + t199 * t234 - t214;
t221 = t229 * t200;
t254 = t193 * t199;
t174 = t191 * t254 - t221;
t188 = 0.1e1 / t192;
t189 = 0.1e1 / t192 ^ 2;
t196 = 0.1e1 / t199 ^ 2;
t231 = t196 * t248;
t195 = 0.1e1 / t199;
t253 = t194 * t195;
t233 = t193 * t253;
t261 = t188 * t195;
t275 = t174 * (t188 * t231 + t189 * t233) + t153 * t261;
t256 = t192 * t199;
t164 = atan2(-t171, t256);
t161 = cos(t164);
t160 = sin(t164);
t268 = t160 * t171;
t152 = t161 * t256 - t268;
t149 = 0.1e1 / t152;
t150 = 0.1e1 / t152 ^ 2;
t274 = 0.2e1 * t174;
t166 = t171 ^ 2;
t260 = t189 * t196;
t165 = t166 * t260 + 0.1e1;
t162 = 0.1e1 / t165;
t212 = t192 * t248 + t194 * t254;
t237 = t171 * t260;
t216 = t200 * t219;
t223 = t229 * t192;
t217 = t194 * t223;
t249 = qJD(5) * t199;
t250 = qJD(1) * t191;
t155 = -t199 * t217 - t191 * t249 - t216 + (t199 * t250 + t200 * t218) * t193;
t239 = t155 * t261;
t141 = (t212 * t237 - t239) * t162;
t208 = -t141 * t171 + t212;
t137 = (-t141 * t256 - t155) * t160 + t208 * t161;
t151 = t149 * t150;
t273 = t137 * t151;
t190 = t188 / t187;
t197 = t195 * t196;
t272 = (t155 * t237 + (-t189 * t197 * t248 - t190 * t196 * t255) * t166) / t165 ^ 2;
t271 = t150 * t174;
t270 = t153 * t150;
t267 = t160 * t174;
t266 = t160 * t192;
t265 = t161 * t171;
t264 = t161 * t174;
t263 = t161 * t193;
t251 = t196 * t200;
t167 = t174 ^ 2;
t147 = t167 * t150 + 0.1e1;
t246 = 0.2e1 * (-t167 * t273 - t174 * t270) / t147 ^ 2;
t245 = 0.2e1 * t276;
t244 = -0.2e1 * t272;
t243 = t151 * t274;
t242 = t188 * t272;
t241 = t150 * t267;
t238 = t171 * t261;
t235 = t189 * t193 * t195;
t232 = t193 * t252;
t228 = t149 * t246;
t227 = t150 * t246;
t226 = t259 * t274;
t225 = t195 * t242;
t211 = t171 * t235 + t229;
t148 = t211 * t162;
t220 = t229 - t148;
t173 = -t191 * t199 + t193 * t221;
t213 = t171 * t251 - t173 * t195;
t210 = t169 * t173 * t191 - t168 * t229;
t157 = 0.1e1 / t159;
t156 = qJD(1) * t175 - t193 * t214 - t200 * t217 - t230;
t145 = 0.1e1 / t147;
t144 = t213 * t188 * t162;
t140 = (-t160 + (t161 * t238 + t160) * t162) * t174;
t139 = -t148 * t265 + (t220 * t266 + t263) * t199;
t138 = t161 * t192 * t200 - t160 * t173 + (-t160 * t256 - t265) * t144;
t136 = t211 * t244 + (t155 * t235 + t250 + (-t188 * t253 + (-t189 * t231 - 0.2e1 * t190 * t233) * t193) * t171) * t162;
t134 = -0.2e1 * t213 * t242 + (-t213 * t189 * t255 + (t155 * t251 - t156 * t195 + (t173 * t251 + (-0.2e1 * t197 * t200 ^ 2 - t195) * t171) * qJD(5)) * t188) * t162;
t133 = (t168 * t191 * t193 + t200 * t236) * t245 + (0.2e1 * t200 * t224 + (t215 + t234) * t168 + ((t154 * t193 + 0.2e1 * t187 * t216) * t191 + (t187 * t249 - 0.2e1 * t192 * t232) * t186) * t169) * t157;
t132 = t139 * t174 * t227 + (-(-t136 * t265 + (t141 * t268 - t155 * t161) * t148) * t271 + (t137 * t243 + t270) * t139 + (-t149 * t259 - (-t148 * t266 + t160 * t223 + t263) * t271) * t248) * t145 + (t228 * t259 + ((-t149 * t258 - (t194 * t220 - t141) * t241) * t193 + (t149 * t219 + (t191 * t137 - (-t136 + t250) * t267 - (t141 * t220 - t194) * t264) * t150) * t192) * t145) * t199;
t1 = [t275 * t162 + t225 * t274, 0, t136, t136, t134, 0; t171 * t228 + (-t155 * t149 + (t137 * t171 + t140 * t153) * t150) * t145 + (t140 * t227 + (0.2e1 * t140 * t273 + (t153 * t162 - t153 - (-t141 * t162 * t238 + t244) * t174) * t150 * t160 + (-(-0.2e1 * t171 * t225 - t141) * t271 + (-(t141 + t239) * t174 + t275 * t171) * t150 * t162) * t161) * t145) * t174, 0, t132, t132 (t138 * t271 - t149 * t175) * t246 + (t138 * t270 + t154 * t149 + (t138 * t243 - t175 * t150) * t137 - (-t192 * t249 + t232 - t134 * t171 - t144 * t155 + (-t144 * t256 - t173) * t141) * t150 * t264 - (-t156 + (-t134 * t199 - t141 * t200) * t192 - t208 * t144) * t241) * t145, 0; t210 * t192 * t245 + (-t210 * t255 + ((qJD(1) * t168 + 0.2e1 * t173 * t269) * t191 + (-t154 * t229 - t156 * t191 + t173 * t219) * t169) * t192) * t157, 0, t133, t133, t169 * t226 * t276 + (t226 * t269 + (t153 * t259 + (-t191 * t255 + t192 * t219) * t174) * t169) * t157, 0;];
JaD_rot  = t1;
