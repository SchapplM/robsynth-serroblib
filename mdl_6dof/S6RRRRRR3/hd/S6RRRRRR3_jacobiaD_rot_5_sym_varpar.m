% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:18
% EndTime: 2019-02-26 22:48:19
% DurationCPUTime: 0.80s
% Computational Cost: add. (4423->98), mult. (4025->206), div. (771->12), fcn. (4686->9), ass. (0->97)
t196 = sin(qJ(1));
t256 = 0.2e1 * t196;
t192 = t196 ^ 2;
t195 = qJ(2) + qJ(3);
t187 = sin(t195);
t182 = t187 ^ 2;
t189 = cos(t195);
t184 = 0.1e1 / t189 ^ 2;
t242 = t182 * t184;
t178 = t192 * t242 + 0.1e1;
t176 = 0.1e1 / t178;
t183 = 0.1e1 / t189;
t197 = cos(qJ(1));
t228 = qJD(1) * t197;
t218 = t187 * t228;
t191 = qJD(2) + qJD(3);
t234 = t191 * t196;
t221 = t184 * t234;
t149 = (-(-t189 * t234 - t218) * t183 + t182 * t221) * t176;
t255 = t149 - t234;
t194 = qJ(4) + qJ(5);
t188 = cos(t194);
t186 = sin(t194);
t232 = t196 * t186;
t235 = t189 * t197;
t171 = t188 * t235 + t232;
t231 = t196 * t187;
t175 = atan2(-t231, -t189);
t173 = cos(t175);
t172 = sin(t175);
t222 = t172 * t231;
t159 = -t173 * t189 - t222;
t156 = 0.1e1 / t159;
t165 = 0.1e1 / t171;
t157 = 0.1e1 / t159 ^ 2;
t166 = 0.1e1 / t171 ^ 2;
t254 = t176 - 0.1e1;
t244 = t173 * t187;
t144 = (-t149 * t196 + t191) * t244 + (t255 * t189 - t218) * t172;
t253 = t144 * t156 * t157;
t190 = qJD(4) + qJD(5);
t207 = t188 * t197 + t189 * t232;
t233 = t191 * t197;
t219 = t187 * t233;
t150 = t207 * qJD(1) - t171 * t190 + t186 * t219;
t230 = t196 * t188;
t170 = t186 * t235 - t230;
t164 = t170 ^ 2;
t163 = t164 * t166 + 0.1e1;
t247 = t166 * t170;
t212 = -qJD(1) * t189 + t190;
t213 = t189 * t190 - qJD(1);
t238 = t186 * t197;
t151 = -t213 * t238 + (t212 * t196 - t219) * t188;
t251 = t151 * t165 * t166;
t252 = (-t150 * t247 - t164 * t251) / t163 ^ 2;
t181 = t187 * t182;
t239 = t183 * t187;
t206 = t191 * (t181 * t183 * t184 + t239);
t240 = t182 * t196;
t210 = t228 * t240;
t250 = (t184 * t210 + t192 * t206) / t178 ^ 2;
t249 = t157 * t187;
t248 = t165 * t186;
t246 = t170 * t188;
t245 = t172 * t196;
t243 = t182 * t183;
t193 = t197 ^ 2;
t241 = t182 * t193;
t237 = t187 * t197;
t236 = t189 * t191;
t229 = qJD(1) * t196;
t154 = t157 * t241 + 0.1e1;
t227 = 0.2e1 * (-t241 * t253 + (t187 * t193 * t236 - t210) * t157) / t154 ^ 2;
t226 = 0.2e1 * t253;
t225 = -0.2e1 * t252;
t224 = t157 * t237;
t223 = t170 * t251;
t217 = 0.1e1 + t242;
t216 = t187 * t227;
t215 = -0.2e1 * t187 * t250;
t214 = t250 * t256;
t211 = t173 * t176 * t243;
t209 = t217 * t197;
t208 = t166 * t246 - t248;
t205 = t191 * t231 + t212 * t197;
t169 = -t189 * t230 + t238;
t162 = t217 * t196 * t176;
t160 = 0.1e1 / t163;
t152 = 0.1e1 / t154;
t148 = (t254 * t187 * t172 - t196 * t211) * t197;
t147 = -t189 * t245 + t244 + (t172 * t189 - t173 * t231) * t162;
t145 = -t217 * t214 + (qJD(1) * t209 + t206 * t256) * t176;
t142 = t225 + 0.2e1 * (-t150 * t160 * t166 + (-t160 * t251 - t166 * t252) * t170) * t170;
t141 = t208 * t225 * t237 + (t208 * t189 * t233 + (-t208 * t229 + ((-t165 * t190 - 0.2e1 * t223) * t188 + (-t150 * t188 + (-t170 * t190 + t151) * t186) * t166) * t197) * t187) * t160;
t140 = (t147 * t249 - t156 * t189) * t197 * t227 + ((-t156 * t229 + (-t147 * t191 - t144) * t197 * t157) * t189 + (-t156 * t233 - (-t145 * t173 * t196 - t255 * t172 + (t149 * t245 - t172 * t191 - t173 * t228) * t162) * t224 + (t157 * t229 + t197 * t226) * t147 - ((t145 - t228) * t172 + ((-t162 * t196 + 0.1e1) * t191 + (t162 - t196) * t149) * t173) * t157 * t235) * t187) * t152;
t1 = [t183 * t197 * t215 + (t191 * t209 - t229 * t239) * t176, t145, t145, 0, 0, 0; (t156 * t216 + (-t156 * t236 + (qJD(1) * t148 + t144) * t249) * t152) * t196 + (t157 * t216 * t148 + (-((t215 - t236 + (t149 * t183 * t240 + t236) * t176) * t172 + (t214 * t243 - t149 * t187 + (-t181 * t221 + (t149 - 0.2e1 * t234) * t187) * t176) * t173) * t224 + (-t157 * t236 + t187 * t226) * t148 + (-t156 + ((-t192 + t193) * t211 + t254 * t222) * t157) * t187 * qJD(1)) * t152) * t197, t140, t140, 0, 0, 0; 0.2e1 * (t165 * t207 + t169 * t247) * t252 + (0.2e1 * t169 * t223 - t213 * t165 * t230 + t205 * t248 + (-t213 * t170 * t232 + t169 * t150 + t151 * t207 - t205 * t246) * t166) * t160, t141, t141, t142, t142, 0;];
JaD_rot  = t1;
