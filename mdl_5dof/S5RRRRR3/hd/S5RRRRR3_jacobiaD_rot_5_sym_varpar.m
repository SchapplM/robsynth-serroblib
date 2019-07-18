% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiaD_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobiaD_rot_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiaD_rot_5_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:38
% EndTime: 2019-07-18 17:19:39
% DurationCPUTime: 0.75s
% Computational Cost: add. (4423->98), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->98)
t197 = sin(qJ(1));
t258 = 0.2e1 * t197;
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
t237 = t192 * t197;
t222 = t185 * t237;
t150 = (-(-t190 * t237 - t219) * t184 + t183 * t222) * t177;
t257 = t150 - t237;
t195 = qJ(4) + qJ(5);
t189 = cos(t195);
t231 = t198 * t189;
t187 = sin(t195);
t235 = t197 * t187;
t172 = t190 * t231 + t235;
t234 = t197 * t188;
t176 = atan2(-t234, -t190);
t174 = cos(t176);
t173 = sin(t176);
t223 = t173 * t234;
t160 = -t174 * t190 - t223;
t157 = 0.1e1 / t160;
t166 = 0.1e1 / t172;
t158 = 0.1e1 / t160 ^ 2;
t167 = 0.1e1 / t172 ^ 2;
t256 = t177 - 0.1e1;
t245 = t174 * t188;
t145 = (-t150 * t197 + t192) * t245 + (t257 * t190 - t219) * t173;
t255 = t145 * t157 * t158;
t191 = qJD(4) + qJD(5);
t208 = t190 * t235 + t231;
t236 = t192 * t198;
t220 = t188 * t236;
t151 = t208 * qJD(1) - t172 * t191 + t187 * t220;
t232 = t198 * t187;
t233 = t197 * t189;
t171 = t190 * t232 - t233;
t165 = t171 ^ 2;
t164 = t165 * t167 + 0.1e1;
t248 = t167 * t171;
t213 = -qJD(1) * t190 + t191;
t214 = t190 * t191 - qJD(1);
t152 = -t214 * t232 + (t213 * t197 - t220) * t189;
t253 = t152 * t166 * t167;
t254 = (-t151 * t248 - t165 * t253) / t164 ^ 2;
t182 = t188 * t183;
t240 = t184 * t188;
t207 = t192 * (t182 * t184 * t185 + t240);
t241 = t183 * t197;
t211 = t229 * t241;
t252 = (t185 * t211 + t193 * t207) / t179 ^ 2;
t251 = t158 * t188;
t250 = t158 * t198;
t249 = t166 * t187;
t247 = t171 * t189;
t246 = t173 * t197;
t244 = t183 * t184;
t194 = t198 ^ 2;
t242 = t183 * t194;
t239 = t188 * t198;
t238 = t190 * t192;
t230 = qJD(1) * t197;
t155 = t158 * t242 + 0.1e1;
t228 = 0.2e1 * (-t242 * t255 + (t188 * t194 * t238 - t211) * t158) / t155 ^ 2;
t227 = 0.2e1 * t255;
t226 = -0.2e1 * t254;
t225 = t158 * t239;
t224 = t171 * t253;
t218 = 0.1e1 + t243;
t217 = t188 * t228;
t216 = -0.2e1 * t188 * t252;
t215 = t252 * t258;
t212 = t174 * t177 * t244;
t210 = t218 * t198;
t209 = t167 * t247 - t249;
t206 = t192 * t234 + t213 * t198;
t170 = -t190 * t233 + t232;
t163 = t218 * t197 * t177;
t161 = 0.1e1 / t164;
t153 = 0.1e1 / t155;
t149 = (t256 * t188 * t173 - t197 * t212) * t198;
t148 = -t190 * t246 + t245 + (t173 * t190 - t174 * t234) * t163;
t146 = -t218 * t215 + (qJD(1) * t210 + t207 * t258) * t177;
t143 = t226 + 0.2e1 * (-t151 * t167 * t161 + (-t161 * t253 - t167 * t254) * t171) * t171;
t142 = t209 * t226 * t239 + (t209 * t190 * t236 + (-t209 * t230 + ((-t166 * t191 - 0.2e1 * t224) * t189 + (-t151 * t189 + (-t171 * t191 + t152) * t187) * t167) * t198) * t188) * t161;
t141 = (t148 * t251 - t157 * t190) * t198 * t228 + ((-t157 * t230 + (-t148 * t192 - t145) * t250) * t190 + (-t157 * t236 - (-t146 * t174 * t197 - t257 * t173 + (t150 * t246 - t173 * t192 - t174 * t229) * t163) * t225 + (t158 * t230 + t198 * t227) * t148 - ((t146 - t229) * t173 + ((-t163 * t197 + 0.1e1) * t192 + (t163 - t197) * t150) * t174) * t190 * t250) * t188) * t153;
t1 = [t198 * t184 * t216 + (t192 * t210 - t230 * t240) * t177, t146, t146, 0, 0; (t157 * t217 + (-t157 * t238 + (qJD(1) * t149 + t145) * t251) * t153) * t197 + (t158 * t217 * t149 + (-((t216 - t238 + (t150 * t184 * t241 + t238) * t177) * t173 + (t215 * t244 - t150 * t188 + (-t182 * t222 + (t150 - 0.2e1 * t237) * t188) * t177) * t174) * t225 + (-t158 * t238 + t188 * t227) * t149 + (-t157 + ((-t193 + t194) * t212 + t256 * t223) * t158) * t188 * qJD(1)) * t153) * t198, t141, t141, 0, 0; 0.2e1 * (t166 * t208 + t170 * t248) * t254 + (0.2e1 * t170 * t224 - t214 * t166 * t233 + t206 * t249 + (-t214 * t171 * t235 + t170 * t151 + t152 * t208 - t206 * t247) * t167) * t161, t142, t142, t143, t143;];
JaD_rot  = t1;
