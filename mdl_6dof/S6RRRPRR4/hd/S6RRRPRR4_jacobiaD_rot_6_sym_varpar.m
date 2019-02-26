% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:42
% EndTime: 2019-02-26 22:17:43
% DurationCPUTime: 0.79s
% Computational Cost: add. (5000->98), mult. (4025->206), div. (771->12), fcn. (4686->9), ass. (0->97)
t196 = sin(qJ(1));
t256 = 0.2e1 * t196;
t193 = t196 ^ 2;
t195 = qJ(2) + qJ(3);
t189 = sin(t195);
t184 = t189 ^ 2;
t190 = cos(t195);
t186 = 0.1e1 / t190 ^ 2;
t241 = t184 * t186;
t178 = t193 * t241 + 0.1e1;
t176 = 0.1e1 / t178;
t185 = 0.1e1 / t190;
t197 = cos(qJ(1));
t228 = qJD(1) * t197;
t218 = t189 * t228;
t192 = qJD(2) + qJD(3);
t234 = t192 * t196;
t220 = t186 * t234;
t151 = (-(-t190 * t234 - t218) * t185 + t184 * t220) * t176;
t255 = t151 - t234;
t188 = pkin(11) + qJ(5) + qJ(6);
t182 = cos(t188);
t181 = sin(t188);
t232 = t196 * t181;
t235 = t190 * t197;
t171 = t182 * t235 + t232;
t230 = t196 * t189;
t175 = atan2(-t230, -t190);
t174 = cos(t175);
t173 = sin(t175);
t222 = t173 * t230;
t162 = -t174 * t190 - t222;
t159 = 0.1e1 / t162;
t165 = 0.1e1 / t171;
t160 = 0.1e1 / t162 ^ 2;
t166 = 0.1e1 / t171 ^ 2;
t254 = t176 - 0.1e1;
t244 = t174 * t189;
t144 = (-t151 * t196 + t192) * t244 + (t255 * t190 - t218) * t173;
t253 = t144 * t159 * t160;
t191 = qJD(5) + qJD(6);
t212 = -qJD(1) * t190 + t191;
t213 = t190 * t191 - qJD(1);
t233 = t192 * t197;
t219 = t189 * t233;
t243 = t181 * t197;
t150 = -t213 * t243 + (t212 * t196 - t219) * t182;
t252 = t150 * t165 * t166;
t207 = t182 * t197 + t190 * t232;
t149 = t207 * qJD(1) - t171 * t191 + t181 * t219;
t231 = t196 * t182;
t170 = t181 * t235 - t231;
t164 = t170 ^ 2;
t157 = t164 * t166 + 0.1e1;
t247 = t166 * t170;
t251 = 0.1e1 / t157 ^ 2 * (-t149 * t247 - t164 * t252);
t183 = t189 * t184;
t238 = t185 * t189;
t206 = t192 * (t183 * t185 * t186 + t238);
t239 = t184 * t196;
t210 = t228 * t239;
t250 = (t186 * t210 + t193 * t206) / t178 ^ 2;
t249 = t160 * t189;
t248 = t165 * t181;
t246 = t170 * t182;
t245 = t173 * t196;
t242 = t184 * t185;
t194 = t197 ^ 2;
t240 = t184 * t194;
t237 = t189 * t197;
t236 = t190 * t192;
t229 = qJD(1) * t196;
t154 = t160 * t240 + 0.1e1;
t227 = 0.2e1 * (-t240 * t253 + (t189 * t194 * t236 - t210) * t160) / t154 ^ 2;
t226 = 0.2e1 * t253;
t225 = -0.2e1 * t251;
t224 = t170 * t252;
t223 = t160 * t237;
t217 = 0.1e1 + t241;
t216 = t189 * t227;
t215 = -0.2e1 * t189 * t250;
t214 = t250 * t256;
t211 = t174 * t176 * t242;
t209 = t217 * t197;
t208 = t166 * t246 - t248;
t205 = t192 * t230 + t212 * t197;
t169 = -t190 * t231 + t243;
t163 = t217 * t196 * t176;
t155 = 0.1e1 / t157;
t152 = 0.1e1 / t154;
t148 = (t254 * t189 * t173 - t196 * t211) * t197;
t147 = -t190 * t245 + t244 + (t173 * t190 - t174 * t230) * t163;
t145 = -t217 * t214 + (qJD(1) * t209 + t206 * t256) * t176;
t142 = t225 + 0.2e1 * (-t149 * t166 * t155 + (-t155 * t252 - t166 * t251) * t170) * t170;
t141 = t208 * t225 * t237 + (t208 * t190 * t233 + (-t208 * t229 + ((-t165 * t191 - 0.2e1 * t224) * t182 + (-t149 * t182 + (-t170 * t191 + t150) * t181) * t166) * t197) * t189) * t155;
t140 = (t147 * t249 - t159 * t190) * t197 * t227 + ((-t159 * t229 + (-t147 * t192 - t144) * t197 * t160) * t190 + (-t159 * t233 - (-t145 * t174 * t196 - t255 * t173 + (t151 * t245 - t173 * t192 - t174 * t228) * t163) * t223 + (t160 * t229 + t197 * t226) * t147 - ((t145 - t228) * t173 + ((-t163 * t196 + 0.1e1) * t192 + (t163 - t196) * t151) * t174) * t160 * t235) * t189) * t152;
t1 = [t197 * t185 * t215 + (t192 * t209 - t229 * t238) * t176, t145, t145, 0, 0, 0; (t159 * t216 + (-t159 * t236 + (qJD(1) * t148 + t144) * t249) * t152) * t196 + (t160 * t216 * t148 + (-((t215 - t236 + (t151 * t185 * t239 + t236) * t176) * t173 + (t214 * t242 - t151 * t189 + (-t183 * t220 + (t151 - 0.2e1 * t234) * t189) * t176) * t174) * t223 + (-t160 * t236 + t189 * t226) * t148 + (-t159 + ((-t193 + t194) * t211 + t254 * t222) * t160) * t189 * qJD(1)) * t152) * t197, t140, t140, 0, 0, 0; 0.2e1 * (t165 * t207 + t169 * t247) * t251 + (0.2e1 * t169 * t224 - t213 * t165 * t231 + t205 * t248 + (-t213 * t170 * t232 + t169 * t149 + t150 * t207 - t205 * t246) * t166) * t155, t141, t141, 0, t142, t142;];
JaD_rot  = t1;
