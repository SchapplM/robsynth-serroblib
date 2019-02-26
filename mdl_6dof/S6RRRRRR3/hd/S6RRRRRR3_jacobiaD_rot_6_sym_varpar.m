% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_rot = S6RRRRRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:18
% EndTime: 2019-02-26 22:48:19
% DurationCPUTime: 0.84s
% Computational Cost: add. (5411->98), mult. (4240->206), div. (789->12), fcn. (4917->9), ass. (0->97)
t199 = sin(qJ(1));
t259 = 0.2e1 * t199;
t196 = t199 ^ 2;
t198 = qJ(2) + qJ(3);
t192 = sin(t198);
t187 = t192 ^ 2;
t193 = cos(t198);
t189 = 0.1e1 / t193 ^ 2;
t244 = t187 * t189;
t182 = t196 * t244 + 0.1e1;
t179 = 0.1e1 / t182;
t188 = 0.1e1 / t193;
t200 = cos(qJ(1));
t231 = qJD(1) * t200;
t221 = t192 * t231;
t195 = qJD(2) + qJD(3);
t237 = t195 * t199;
t223 = t189 * t237;
t154 = (-(-t193 * t237 - t221) * t188 + t187 * t223) * t179;
t258 = t154 - t237;
t194 = qJ(4) + qJ(5) + qJ(6);
t185 = cos(t194);
t184 = sin(t194);
t235 = t199 * t184;
t238 = t193 * t200;
t174 = t185 * t238 + t235;
t233 = t199 * t192;
t178 = atan2(-t233, -t193);
t177 = cos(t178);
t176 = sin(t178);
t225 = t176 * t233;
t165 = -t177 * t193 - t225;
t162 = 0.1e1 / t165;
t168 = 0.1e1 / t174;
t163 = 0.1e1 / t165 ^ 2;
t169 = 0.1e1 / t174 ^ 2;
t257 = t179 - 0.1e1;
t247 = t177 * t192;
t147 = (-t154 * t199 + t195) * t247 + (t258 * t193 - t221) * t176;
t256 = t147 * t162 * t163;
t191 = qJD(4) + qJD(5) + qJD(6);
t210 = t185 * t200 + t193 * t235;
t236 = t195 * t200;
t222 = t192 * t236;
t152 = t210 * qJD(1) - t174 * t191 + t184 * t222;
t234 = t199 * t185;
t173 = t184 * t238 - t234;
t167 = t173 ^ 2;
t161 = t167 * t169 + 0.1e1;
t250 = t169 * t173;
t215 = -qJD(1) * t193 + t191;
t216 = t191 * t193 - qJD(1);
t246 = t184 * t200;
t153 = -t216 * t246 + (t215 * t199 - t222) * t185;
t254 = t153 * t168 * t169;
t255 = (-t152 * t250 - t167 * t254) / t161 ^ 2;
t186 = t192 * t187;
t241 = t188 * t192;
t209 = t195 * (t186 * t188 * t189 + t241);
t242 = t187 * t199;
t213 = t231 * t242;
t253 = (t189 * t213 + t196 * t209) / t182 ^ 2;
t252 = t163 * t192;
t251 = t168 * t184;
t249 = t173 * t185;
t248 = t176 * t199;
t245 = t187 * t188;
t197 = t200 ^ 2;
t243 = t187 * t197;
t240 = t192 * t200;
t239 = t193 * t195;
t232 = qJD(1) * t199;
t157 = t163 * t243 + 0.1e1;
t230 = 0.2e1 * (-t243 * t256 + (t192 * t197 * t239 - t213) * t163) / t157 ^ 2;
t229 = 0.2e1 * t256;
t228 = -0.2e1 * t255;
t227 = t173 * t254;
t226 = t163 * t240;
t220 = 0.1e1 + t244;
t219 = t192 * t230;
t218 = -0.2e1 * t192 * t253;
t217 = t253 * t259;
t214 = t177 * t179 * t245;
t212 = t220 * t200;
t211 = t169 * t249 - t251;
t208 = t195 * t233 + t215 * t200;
t172 = -t193 * t234 + t246;
t166 = t220 * t199 * t179;
t158 = 0.1e1 / t161;
t155 = 0.1e1 / t157;
t151 = (t257 * t192 * t176 - t199 * t214) * t200;
t150 = -t193 * t248 + t247 + (t176 * t193 - t177 * t233) * t166;
t148 = -t220 * t217 + (qJD(1) * t212 + t209 * t259) * t179;
t145 = t228 + 0.2e1 * (-t152 * t158 * t169 + (-t158 * t254 - t169 * t255) * t173) * t173;
t144 = t211 * t228 * t240 + (t211 * t193 * t236 + (-t211 * t232 + ((-t168 * t191 - 0.2e1 * t227) * t185 + (-t152 * t185 + (-t173 * t191 + t153) * t184) * t169) * t200) * t192) * t158;
t143 = (t150 * t252 - t162 * t193) * t200 * t230 + ((-t162 * t232 + (-t150 * t195 - t147) * t200 * t163) * t193 + (-t162 * t236 - (-t148 * t177 * t199 - t258 * t176 + (t154 * t248 - t176 * t195 - t177 * t231) * t166) * t226 + (t163 * t232 + t200 * t229) * t150 - ((t148 - t231) * t176 + ((-t166 * t199 + 0.1e1) * t195 + (t166 - t199) * t154) * t177) * t163 * t238) * t192) * t155;
t1 = [t188 * t200 * t218 + (t195 * t212 - t232 * t241) * t179, t148, t148, 0, 0, 0; (t162 * t219 + (-t162 * t239 + (qJD(1) * t151 + t147) * t252) * t155) * t199 + (t163 * t219 * t151 + (-((t218 - t239 + (t154 * t188 * t242 + t239) * t179) * t176 + (t217 * t245 - t154 * t192 + (-t186 * t223 + (t154 - 0.2e1 * t237) * t192) * t179) * t177) * t226 + (-t163 * t239 + t192 * t229) * t151 + (-t162 + ((-t196 + t197) * t214 + t257 * t225) * t163) * t192 * qJD(1)) * t155) * t200, t143, t143, 0, 0, 0; 0.2e1 * (t168 * t210 + t172 * t250) * t255 + (0.2e1 * t172 * t227 - t216 * t168 * t234 + t208 * t251 + (-t216 * t173 * t235 + t172 * t152 + t153 * t210 - t208 * t249) * t169) * t158, t144, t144, t145, t145, t145;];
JaD_rot  = t1;
