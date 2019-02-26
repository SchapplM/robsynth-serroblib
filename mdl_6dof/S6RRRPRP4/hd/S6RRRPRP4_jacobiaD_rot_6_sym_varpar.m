% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:05
% EndTime: 2019-02-26 22:11:06
% DurationCPUTime: 1.03s
% Computational Cost: add. (5529->122), mult. (8382->269), div. (1515->15), fcn. (10508->9), ass. (0->120)
t194 = qJ(2) + qJ(3);
t187 = sin(t194);
t196 = sin(qJ(1));
t197 = cos(qJ(5));
t244 = t196 * t197;
t195 = sin(qJ(5));
t198 = cos(qJ(1));
t245 = t195 * t198;
t174 = t187 * t245 + t244;
t170 = 0.1e1 / t174 ^ 2;
t188 = cos(t194);
t183 = t188 ^ 2;
t193 = t198 ^ 2;
t256 = t183 * t193;
t229 = t170 * t256;
t166 = 0.1e1 + t229;
t216 = qJD(5) * t187 + qJD(1);
t212 = t216 * t197;
t215 = qJD(1) * t187 + qJD(5);
t189 = qJD(2) + qJD(3);
t248 = t189 * t198;
t222 = t188 * t248;
t158 = t198 * t212 + (-t196 * t215 + t222) * t195;
t169 = 0.1e1 / t174;
t263 = t158 * t169 * t170;
t214 = t256 * t263;
t239 = qJD(1) * t198;
t220 = t196 * t239;
t224 = t188 * t189 * t193;
t272 = (-t214 + (-t183 * t220 - t187 * t224) * t170) / t166 ^ 2;
t184 = 0.1e1 / t188;
t175 = t187 * t244 + t245;
t243 = t197 * t198;
t246 = t195 * t196;
t176 = -t187 * t246 + t243;
t190 = 0.1e1 / t197;
t191 = 0.1e1 / t197 ^ 2;
t247 = t191 * t195;
t209 = t175 * t247 + t176 * t190;
t271 = t184 * t209;
t250 = t188 * t198;
t251 = t188 * t197;
t165 = atan2(t175, t251);
t160 = cos(t165);
t159 = sin(t165);
t261 = t159 * t175;
t154 = t160 * t251 + t261;
t151 = 0.1e1 / t154;
t152 = 0.1e1 / t154 ^ 2;
t221 = t187 * t243;
t173 = -t221 + t246;
t270 = 0.2e1 * t173;
t269 = -0.2e1 * t195;
t168 = t173 ^ 2;
t149 = t152 * t168 + 0.1e1;
t157 = qJD(1) * t175 + qJD(5) * t174 - t197 * t222;
t264 = t157 * t152;
t172 = t175 ^ 2;
t185 = 0.1e1 / t188 ^ 2;
t254 = t185 * t191;
t167 = t172 * t254 + 0.1e1;
t163 = 0.1e1 / t167;
t238 = qJD(5) * t195;
t253 = t187 * t189;
t208 = -t188 * t238 - t197 * t253;
t227 = t175 * t254;
t252 = t188 * t196;
t223 = t189 * t252;
t237 = qJD(5) * t197;
t155 = -qJD(1) * t221 - t197 * t223 - t198 * t237 + t216 * t246;
t255 = t184 * t190;
t230 = t155 * t255;
t143 = (-t208 * t227 - t230) * t163;
t207 = -t143 * t175 - t208;
t139 = (-t143 * t251 - t155) * t159 - t207 * t160;
t153 = t151 * t152;
t267 = t139 * t153;
t268 = (-t168 * t267 + t173 * t264) / t149 ^ 2;
t186 = t184 / t183;
t192 = t190 * t191;
t266 = (-t155 * t227 + (t185 * t192 * t238 + t186 * t191 * t253) * t172) / t167 ^ 2;
t265 = t152 * t173;
t262 = t159 * t173;
t260 = t159 * t188;
t259 = t160 * t173;
t258 = t160 * t175;
t257 = t160 * t187;
t249 = t189 * t190;
t242 = t198 * t151;
t226 = t185 * t187 * t190;
t211 = t175 * t226 + t196;
t150 = t211 * t163;
t241 = t150 - t196;
t240 = qJD(1) * t196;
t236 = 0.2e1 * t268;
t235 = 0.2e1 * t272;
t234 = -0.2e1 * t266;
t233 = t153 * t270;
t232 = t151 * t268;
t231 = t152 * t262;
t228 = t175 * t255;
t225 = t187 * t249;
t219 = t191 * t238;
t218 = t152 * t236;
t217 = t250 * t270;
t213 = 0.2e1 * t255 * t266;
t210 = -t170 * t176 * t198 - t169 * t196;
t206 = t157 * t255 - (-t184 * t219 - t185 * t225) * t173;
t161 = 0.1e1 / t166;
t156 = t196 * t212 + (t198 * t215 + t223) * t195;
t147 = 0.1e1 / t149;
t146 = t163 * t271;
t142 = (-t159 + (-t160 * t228 + t159) * t163) * t173;
t141 = t150 * t258 + (-t241 * t260 - t257) * t197;
t140 = -t160 * t188 * t195 + t159 * t176 + (-t159 * t251 + t258) * t146;
t138 = t211 * t234 + (-t155 * t226 + t239 + (t184 * t249 + (t185 * t219 + 0.2e1 * t186 * t225) * t187) * t175) * t163;
t136 = t234 * t271 + (t209 * t185 * t253 + (-t155 * t247 - t156 * t190 + (t176 * t247 + (0.2e1 * t192 * t195 ^ 2 + t190) * t175) * qJD(5)) * t184) * t163;
t135 = (-t169 * t187 * t198 - t195 * t229) * t235 + (t214 * t269 + (-t187 * t240 + t222) * t169 + ((-t158 * t198 + t224 * t269) * t187 + (t193 * t237 + t220 * t269) * t183) * t170) * t161;
t134 = t141 * t173 * t218 + (-(t138 * t258 + (-t143 * t261 - t155 * t160) * t150) * t265 + (t139 * t233 - t264) * t141 + (t188 * t242 - (t150 * t260 - t159 * t252 + t257) * t265) * t238) * t147 + (0.2e1 * t232 * t250 + ((t189 * t242 - (t189 * t241 + t143) * t231) * t187 + (t151 * t240 + (t198 * t139 - (-t138 + t239) * t262 - (-t143 * t241 - t189) * t259) * t152) * t188) * t147) * t197;
t1 = [-t163 * t206 + t173 * t213, t138, t138, 0, t136, 0; -0.2e1 * t175 * t232 + (-t155 * t151 + (-t139 * t175 - t142 * t157) * t152) * t147 + (t142 * t218 + (0.2e1 * t142 * t267 + (-t157 * t163 + t157 - (t143 * t163 * t228 + t234) * t173) * t152 * t159 + (-(t175 * t213 - t143) * t265 + (-(t143 + t230) * t173 + t206 * t175) * t152 * t163) * t160) * t147) * t173, t134, t134, 0 (t140 * t265 - t151 * t174) * t236 + (-t140 * t264 + t158 * t151 + (t140 * t233 - t174 * t152) * t139 - (-t188 * t237 + t195 * t253 + t136 * t175 - t146 * t155 + (-t146 * t251 + t176) * t143) * t152 * t259 - (-t156 + (-t136 * t197 + t143 * t195) * t188 + t207 * t146) * t231) * t147, 0; t210 * t188 * t235 + (t210 * t253 + ((qJD(1) * t169 - 0.2e1 * t176 * t263) * t198 + (-t156 * t198 + (-qJD(1) * t176 - t158) * t196) * t170) * t188) * t161, t135, t135, 0, t170 * t217 * t272 + (t217 * t263 + (-t157 * t250 + (t187 * t248 + t188 * t240) * t173) * t170) * t161, 0;];
JaD_rot  = t1;
