% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:17
% EndTime: 2019-02-26 22:32:18
% DurationCPUTime: 0.72s
% Computational Cost: add. (5000->98), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->98)
t196 = sin(qJ(1));
t257 = 0.2e1 * t196;
t193 = t196 ^ 2;
t195 = qJ(2) + qJ(3);
t189 = sin(t195);
t184 = t189 ^ 2;
t190 = cos(t195);
t186 = 0.1e1 / t190 ^ 2;
t242 = t184 * t186;
t178 = t193 * t242 + 0.1e1;
t176 = 0.1e1 / t178;
t185 = 0.1e1 / t190;
t197 = cos(qJ(1));
t228 = qJD(1) * t197;
t218 = t189 * t228;
t192 = qJD(2) + qJD(3);
t236 = t192 * t196;
t221 = t186 * t236;
t151 = (-(-t190 * t236 - t218) * t185 + t184 * t221) * t176;
t256 = t151 - t236;
t188 = qJ(4) + pkin(11) + qJ(6);
t182 = cos(t188);
t230 = t197 * t182;
t181 = sin(t188);
t234 = t196 * t181;
t171 = t190 * t230 + t234;
t232 = t196 * t189;
t175 = atan2(-t232, -t190);
t174 = cos(t175);
t173 = sin(t175);
t222 = t173 * t232;
t162 = -t174 * t190 - t222;
t159 = 0.1e1 / t162;
t165 = 0.1e1 / t171;
t160 = 0.1e1 / t162 ^ 2;
t166 = 0.1e1 / t171 ^ 2;
t255 = t176 - 0.1e1;
t244 = t174 * t189;
t144 = (-t151 * t196 + t192) * t244 + (t256 * t190 - t218) * t173;
t254 = t144 * t159 * t160;
t191 = qJD(4) + qJD(6);
t207 = t190 * t234 + t230;
t235 = t192 * t197;
t219 = t189 * t235;
t149 = t207 * qJD(1) - t171 * t191 + t181 * t219;
t231 = t197 * t181;
t233 = t196 * t182;
t170 = t190 * t231 - t233;
t164 = t170 ^ 2;
t157 = t164 * t166 + 0.1e1;
t247 = t166 * t170;
t212 = -qJD(1) * t190 + t191;
t213 = t190 * t191 - qJD(1);
t150 = -t213 * t231 + (t212 * t196 - t219) * t182;
t252 = t150 * t165 * t166;
t253 = (-t149 * t247 - t164 * t252) / t157 ^ 2;
t183 = t189 * t184;
t239 = t185 * t189;
t206 = t192 * (t183 * t185 * t186 + t239);
t240 = t184 * t196;
t210 = t228 * t240;
t251 = (t186 * t210 + t193 * t206) / t178 ^ 2;
t250 = t160 * t189;
t249 = t160 * t197;
t248 = t165 * t181;
t246 = t170 * t182;
t245 = t173 * t196;
t243 = t184 * t185;
t194 = t197 ^ 2;
t241 = t184 * t194;
t238 = t189 * t197;
t237 = t190 * t192;
t229 = qJD(1) * t196;
t154 = t160 * t241 + 0.1e1;
t227 = 0.2e1 * (-t241 * t254 + (t189 * t194 * t237 - t210) * t160) / t154 ^ 2;
t226 = 0.2e1 * t254;
t225 = -0.2e1 * t253;
t224 = t160 * t238;
t223 = t170 * t252;
t217 = 0.1e1 + t242;
t216 = t189 * t227;
t215 = -0.2e1 * t189 * t251;
t214 = t251 * t257;
t211 = t174 * t176 * t243;
t209 = t217 * t197;
t208 = t166 * t246 - t248;
t205 = t192 * t232 + t212 * t197;
t169 = -t190 * t233 + t231;
t163 = t217 * t196 * t176;
t155 = 0.1e1 / t157;
t152 = 0.1e1 / t154;
t148 = (t255 * t189 * t173 - t196 * t211) * t197;
t147 = -t190 * t245 + t244 + (t173 * t190 - t174 * t232) * t163;
t145 = -t217 * t214 + (qJD(1) * t209 + t206 * t257) * t176;
t142 = t225 + 0.2e1 * (-t149 * t166 * t155 + (-t155 * t252 - t166 * t253) * t170) * t170;
t141 = t208 * t225 * t238 + (t208 * t190 * t235 + (-t208 * t229 + ((-t165 * t191 - 0.2e1 * t223) * t182 + (-t149 * t182 + (-t170 * t191 + t150) * t181) * t166) * t197) * t189) * t155;
t140 = (t147 * t250 - t159 * t190) * t197 * t227 + ((-t159 * t229 + (-t147 * t192 - t144) * t249) * t190 + (-t159 * t235 - (-t145 * t174 * t196 - t256 * t173 + (t151 * t245 - t173 * t192 - t174 * t228) * t163) * t224 + (t160 * t229 + t197 * t226) * t147 - ((t145 - t228) * t173 + ((-t163 * t196 + 0.1e1) * t192 + (t163 - t196) * t151) * t174) * t190 * t249) * t189) * t152;
t1 = [t197 * t185 * t215 + (t192 * t209 - t229 * t239) * t176, t145, t145, 0, 0, 0; (t159 * t216 + (-t159 * t237 + (qJD(1) * t148 + t144) * t250) * t152) * t196 + (t160 * t216 * t148 + (-((t215 - t237 + (t151 * t185 * t240 + t237) * t176) * t173 + (t214 * t243 - t151 * t189 + (-t183 * t221 + (t151 - 0.2e1 * t236) * t189) * t176) * t174) * t224 + (-t160 * t237 + t189 * t226) * t148 + (-t159 + ((-t193 + t194) * t211 + t255 * t222) * t160) * t189 * qJD(1)) * t152) * t197, t140, t140, 0, 0, 0; 0.2e1 * (t165 * t207 + t169 * t247) * t253 + (0.2e1 * t169 * t223 - t213 * t165 * t233 + t205 * t248 + (-t213 * t170 * t234 + t169 * t149 + t150 * t207 - t205 * t246) * t166) * t155, t141, t141, t142, 0, t142;];
JaD_rot  = t1;
