% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP8
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
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:09
% EndTime: 2019-02-26 22:13:10
% DurationCPUTime: 0.98s
% Computational Cost: add. (1265->118), mult. (4276->252), div. (493->12), fcn. (5073->11), ass. (0->111)
t191 = sin(qJ(2));
t183 = t191 ^ 2;
t195 = cos(qJ(2));
t186 = 0.1e1 / t195 ^ 2;
t248 = t183 * t186;
t192 = sin(qJ(1));
t184 = t192 ^ 2;
t178 = t184 * t248 + 0.1e1;
t185 = 0.1e1 / t195;
t245 = t185 * t191;
t268 = t191 * t248;
t204 = qJD(2) * (t185 * t268 + t245);
t196 = cos(qJ(1));
t237 = qJD(1) * t196;
t246 = t183 * t192;
t213 = t237 * t246;
t251 = (t184 * t204 + t186 * t213) / t178 ^ 2;
t269 = -0.2e1 * t251;
t220 = 0.1e1 + t248;
t267 = t192 * t220;
t194 = cos(qJ(3));
t266 = (-qJD(3) + qJD(5)) * t194;
t190 = sin(qJ(3));
t233 = qJD(3) * t190;
t265 = -qJD(5) * t190 + t233;
t239 = t196 * t194;
t241 = t192 * t195;
t205 = t190 * t241 + t239;
t234 = qJD(2) * t196;
t221 = t191 * t234;
t223 = t195 * t239;
t143 = t205 * qJD(1) - qJD(3) * t223 + t190 * t221 - t192 * t233;
t215 = -qJD(1) * t195 + qJD(3);
t216 = qJD(3) * t195 - qJD(1);
t240 = t196 * t190;
t144 = -t216 * t240 + (t215 * t192 - t221) * t194;
t189 = sin(qJ(5));
t193 = cos(qJ(5));
t242 = t192 * t194;
t171 = t195 * t240 - t242;
t172 = t192 * t190 + t223;
t209 = t171 * t193 - t172 * t189;
t136 = t209 * qJD(5) - t143 * t189 + t144 * t193;
t159 = t171 * t189 + t172 * t193;
t151 = 0.1e1 / t159;
t207 = t189 * t190 + t193 * t194;
t208 = t189 * t194 - t190 * t193;
t152 = 0.1e1 / t159 ^ 2;
t255 = t152 * t209;
t264 = t208 * t151 + t207 * t255;
t243 = t192 * t191;
t177 = atan2(t243, t195);
t174 = cos(t177);
t173 = sin(t177);
t225 = t173 * t243;
t163 = t174 * t195 + t225;
t160 = 0.1e1 / t163;
t161 = 0.1e1 / t163 ^ 2;
t263 = 0.2e1 * t191;
t175 = 0.1e1 / t178;
t262 = t175 - 0.1e1;
t135 = t159 * qJD(5) + t143 * t193 + t144 * t189;
t150 = t209 ^ 2;
t141 = t150 * t152 + 0.1e1;
t153 = t151 * t152;
t258 = t136 * t153;
t261 = (-t135 * t255 - t150 * t258) / t141 ^ 2;
t188 = t196 ^ 2;
t247 = t183 * t188;
t149 = t161 * t247 + 0.1e1;
t235 = qJD(2) * t195;
t222 = t191 * t237;
t236 = qJD(2) * t192;
t142 = ((t192 * t235 + t222) * t185 + t236 * t248) * t175;
t249 = t174 * t191;
t133 = (t142 * t192 - qJD(2)) * t249 + (t222 + (-t142 + t236) * t195) * t173;
t259 = t133 * t160 * t161;
t260 = (-t247 * t259 + (t188 * t191 * t235 - t213) * t161) / t149 ^ 2;
t257 = t142 * t173;
t256 = t142 * t191;
t244 = t191 * t196;
t167 = t207 * t244;
t254 = t152 * t167;
t253 = t161 * t191;
t252 = t161 * t196;
t165 = t175 * t267;
t250 = t165 * t192;
t238 = qJD(1) * t192;
t229 = 0.2e1 * t261;
t228 = -0.2e1 * t259;
t227 = -0.2e1 * t153 * t209;
t226 = t161 * t244;
t224 = t175 * t183 * t185;
t219 = -0.2e1 * t191 * t260;
t218 = t136 * t227;
t217 = t185 * t269;
t214 = t192 * t224;
t212 = t220 * t196;
t170 = -t194 * t241 + t240;
t210 = -t170 * t189 - t193 * t205;
t155 = t170 * t193 - t189 * t205;
t206 = t215 * t196;
t166 = t208 * t244;
t147 = 0.1e1 / t149;
t146 = t194 * t206 + (qJD(2) * t191 * t194 + t216 * t190) * t192;
t145 = -t216 * t242 + (t191 * t236 + t206) * t190;
t139 = 0.1e1 / t141;
t138 = (-t262 * t191 * t173 + t174 * t214) * t196;
t137 = t173 * t241 - t249 + (-t173 * t195 + t174 * t243) * t165;
t134 = t267 * t269 + (qJD(1) * t212 + 0.2e1 * t192 * t204) * t175;
t1 = [t217 * t244 + (qJD(2) * t212 - t238 * t245) * t175, t134, 0, 0, 0, 0; (t160 * t219 + (t160 * t235 + (-qJD(1) * t138 - t133) * t253) * t147) * t192 + (t161 * t219 * t138 + (((-t142 * t214 - t262 * t235 + t251 * t263) * t173 + (t217 * t246 + t256 + (-t256 + (t263 + t268) * t236) * t175) * t174) * t226 + (t161 * t235 + t191 * t228) * t138 + (t160 + ((-t184 + t188) * t174 * t224 + t262 * t225) * t161) * t191 * qJD(1)) * t147) * t196, 0.2e1 * (-t137 * t253 + t160 * t195) * t196 * t260 + ((t160 * t238 + (qJD(2) * t137 + t133) * t252) * t195 + (t160 * t234 + (t134 * t174 * t192 - t173 * t236 - t250 * t257 + t257 + (qJD(2) * t173 + t174 * t237) * t165) * t226 + (-t161 * t238 + t196 * t228) * t137 + ((-t134 + t237) * t173 + ((-0.1e1 + t250) * qJD(2) + (-t165 + t192) * t142) * t174) * t195 * t252) * t191) * t147, 0, 0, 0, 0; (t151 * t210 - t155 * t255) * t229 + ((t155 * qJD(5) - t145 * t193 + t146 * t189) * t151 + t155 * t218 + (t210 * t136 + (t210 * qJD(5) + t145 * t189 + t146 * t193) * t209 - t155 * t135) * t152) * t139 (t151 * t166 + t209 * t254) * t229 + (t135 * t254 + (t166 * t152 - t167 * t227) * t136 - t264 * t195 * t234 + (t264 * t238 + ((-t266 * t151 + t265 * t255) * t193 + (t265 * t151 + t266 * t255) * t189) * t196) * t191) * t139 (t151 * t159 + t209 * t255) * t229 + (-t136 * t151 - t209 * t218 + (0.2e1 * t209 * t135 + t159 * t136) * t152) * t139, 0, -0.2e1 * t261 - 0.2e1 * (t135 * t152 * t139 - (-t139 * t258 - t152 * t261) * t209) * t209, 0;];
JaD_rot  = t1;
