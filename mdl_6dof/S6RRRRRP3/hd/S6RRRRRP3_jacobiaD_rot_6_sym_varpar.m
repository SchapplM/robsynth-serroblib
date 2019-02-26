% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:42
% EndTime: 2019-02-26 22:40:43
% DurationCPUTime: 0.82s
% Computational Cost: add. (4423->98), mult. (4025->206), div. (771->12), fcn. (4686->9), ass. (0->97)
t197 = sin(qJ(1));
t257 = 0.2e1 * t197;
t193 = t197 ^ 2;
t196 = qJ(2) + qJ(3);
t188 = sin(t196);
t183 = t188 ^ 2;
t190 = cos(t196);
t185 = 0.1e1 / t190 ^ 2;
t243 = t183 * t185;
t179 = t193 * t243 + 0.1e1;
t177 = 0.1e1 / t179;
t184 = 0.1e1 / t190;
t198 = cos(qJ(1));
t229 = qJD(1) * t198;
t219 = t188 * t229;
t192 = qJD(2) + qJD(3);
t235 = t192 * t197;
t222 = t185 * t235;
t150 = (-(-t190 * t235 - t219) * t184 + t183 * t222) * t177;
t256 = t150 - t235;
t195 = qJ(4) + qJ(5);
t189 = cos(t195);
t187 = sin(t195);
t233 = t197 * t187;
t236 = t190 * t198;
t172 = t189 * t236 + t233;
t232 = t197 * t188;
t176 = atan2(-t232, -t190);
t174 = cos(t176);
t173 = sin(t176);
t223 = t173 * t232;
t160 = -t174 * t190 - t223;
t157 = 0.1e1 / t160;
t166 = 0.1e1 / t172;
t158 = 0.1e1 / t160 ^ 2;
t167 = 0.1e1 / t172 ^ 2;
t255 = t177 - 0.1e1;
t245 = t174 * t188;
t145 = (-t150 * t197 + t192) * t245 + (t256 * t190 - t219) * t173;
t254 = t145 * t157 * t158;
t191 = qJD(4) + qJD(5);
t208 = t189 * t198 + t190 * t233;
t234 = t192 * t198;
t221 = t188 * t234;
t151 = t208 * qJD(1) - t172 * t191 + t187 * t221;
t231 = t197 * t189;
t171 = t187 * t236 - t231;
t165 = t171 ^ 2;
t164 = t165 * t167 + 0.1e1;
t248 = t167 * t171;
t213 = -qJD(1) * t190 + t191;
t214 = t190 * t191 - qJD(1);
t239 = t187 * t198;
t152 = -t214 * t239 + (t213 * t197 - t221) * t189;
t252 = t152 * t166 * t167;
t253 = (-t151 * t248 - t165 * t252) / t164 ^ 2;
t182 = t188 * t183;
t240 = t184 * t188;
t207 = t192 * (t182 * t184 * t185 + t240);
t241 = t183 * t197;
t211 = t229 * t241;
t251 = (t185 * t211 + t193 * t207) / t179 ^ 2;
t250 = t158 * t188;
t249 = t166 * t187;
t247 = t171 * t189;
t246 = t173 * t197;
t244 = t183 * t184;
t194 = t198 ^ 2;
t242 = t183 * t194;
t238 = t188 * t198;
t237 = t190 * t192;
t230 = qJD(1) * t197;
t155 = t158 * t242 + 0.1e1;
t228 = 0.2e1 * (-t242 * t254 + (t188 * t194 * t237 - t211) * t158) / t155 ^ 2;
t227 = 0.2e1 * t254;
t226 = -0.2e1 * t253;
t225 = t171 * t252;
t224 = t158 * t238;
t218 = 0.1e1 + t243;
t217 = t188 * t228;
t216 = -0.2e1 * t188 * t251;
t215 = t251 * t257;
t212 = t174 * t177 * t244;
t210 = t218 * t198;
t209 = t167 * t247 - t249;
t206 = t192 * t232 + t213 * t198;
t170 = -t190 * t231 + t239;
t163 = t218 * t197 * t177;
t161 = 0.1e1 / t164;
t153 = 0.1e1 / t155;
t149 = (t255 * t188 * t173 - t197 * t212) * t198;
t148 = -t190 * t246 + t245 + (t173 * t190 - t174 * t232) * t163;
t146 = -t218 * t215 + (qJD(1) * t210 + t207 * t257) * t177;
t143 = t226 + 0.2e1 * (-t151 * t161 * t167 + (-t161 * t252 - t167 * t253) * t171) * t171;
t142 = t209 * t226 * t238 + (t209 * t190 * t234 + (-t209 * t230 + ((-t166 * t191 - 0.2e1 * t225) * t189 + (-t151 * t189 + (-t171 * t191 + t152) * t187) * t167) * t198) * t188) * t161;
t141 = (t148 * t250 - t157 * t190) * t198 * t228 + ((-t157 * t230 + (-t148 * t192 - t145) * t198 * t158) * t190 + (-t157 * t234 - (-t146 * t174 * t197 - t256 * t173 + (t150 * t246 - t173 * t192 - t174 * t229) * t163) * t224 + (t158 * t230 + t198 * t227) * t148 - ((t146 - t229) * t173 + ((-t163 * t197 + 0.1e1) * t192 + (t163 - t197) * t150) * t174) * t158 * t236) * t188) * t153;
t1 = [t184 * t198 * t216 + (t192 * t210 - t230 * t240) * t177, t146, t146, 0, 0, 0; (t157 * t217 + (-t157 * t237 + (qJD(1) * t149 + t145) * t250) * t153) * t197 + (t158 * t217 * t149 + (-((t216 - t237 + (t150 * t184 * t241 + t237) * t177) * t173 + (t215 * t244 - t150 * t188 + (-t182 * t222 + (t150 - 0.2e1 * t235) * t188) * t177) * t174) * t224 + (-t158 * t237 + t188 * t227) * t149 + (-t157 + ((-t193 + t194) * t212 + t255 * t223) * t158) * t188 * qJD(1)) * t153) * t198, t141, t141, 0, 0, 0; 0.2e1 * (t166 * t208 + t170 * t248) * t253 + (0.2e1 * t170 * t225 - t214 * t166 * t231 + t206 * t249 + (-t214 * t171 * t233 + t170 * t151 + t152 * t208 - t206 * t247) * t167) * t161, t142, t142, t143, t143, 0;];
JaD_rot  = t1;
