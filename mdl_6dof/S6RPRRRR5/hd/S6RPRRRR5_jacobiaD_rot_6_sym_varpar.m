% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:10
% EndTime: 2019-02-26 21:17:11
% DurationCPUTime: 0.71s
% Computational Cost: add. (6414->98), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->98)
t197 = sin(qJ(1));
t258 = 0.2e1 * t197;
t194 = t197 ^ 2;
t189 = pkin(11) + qJ(3) + qJ(4);
t187 = sin(t189);
t183 = t187 ^ 2;
t188 = cos(t189);
t185 = 0.1e1 / t188 ^ 2;
t243 = t183 * t185;
t179 = t194 * t243 + 0.1e1;
t176 = 0.1e1 / t179;
t184 = 0.1e1 / t188;
t198 = cos(qJ(1));
t229 = qJD(1) * t198;
t219 = t187 * t229;
t193 = qJD(3) + qJD(4);
t235 = t193 * t197;
t222 = t185 * t235;
t150 = (-(-t188 * t235 - t219) * t184 + t183 * t222) * t176;
t257 = t150 - t235;
t196 = qJ(5) + qJ(6);
t190 = sin(t196);
t232 = t197 * t190;
t191 = cos(t196);
t236 = t191 * t198;
t172 = t188 * t236 + t232;
t233 = t197 * t187;
t175 = atan2(-t233, -t188);
t174 = cos(t175);
t173 = sin(t175);
t223 = t173 * t233;
t160 = -t174 * t188 - t223;
t157 = 0.1e1 / t160;
t166 = 0.1e1 / t172;
t158 = 0.1e1 / t160 ^ 2;
t167 = 0.1e1 / t172 ^ 2;
t256 = t176 - 0.1e1;
t245 = t174 * t187;
t145 = (-t150 * t197 + t193) * t245 + (t257 * t188 - t219) * t173;
t255 = t145 * t157 * t158;
t192 = qJD(5) + qJD(6);
t213 = -qJD(1) * t188 + t192;
t214 = t188 * t192 - qJD(1);
t234 = t193 * t198;
t220 = t187 * t234;
t237 = t190 * t198;
t152 = -t214 * t237 + (t213 * t197 - t220) * t191;
t254 = t152 * t166 * t167;
t182 = t187 * t183;
t240 = t184 * t187;
t207 = t193 * (t182 * t184 * t185 + t240);
t241 = t183 * t197;
t211 = t229 * t241;
t253 = (t185 * t211 + t194 * t207) / t179 ^ 2;
t252 = t158 * t187;
t251 = t158 * t198;
t208 = t188 * t232 + t236;
t151 = t208 * qJD(1) - t172 * t192 + t190 * t220;
t231 = t197 * t191;
t171 = t188 * t237 - t231;
t165 = t171 ^ 2;
t164 = t165 * t167 + 0.1e1;
t248 = t167 * t171;
t250 = 0.1e1 / t164 ^ 2 * (-t151 * t248 - t165 * t254);
t249 = t166 * t190;
t247 = t171 * t191;
t246 = t173 * t197;
t244 = t183 * t184;
t195 = t198 ^ 2;
t242 = t183 * t195;
t239 = t187 * t198;
t238 = t188 * t193;
t230 = qJD(1) * t197;
t155 = t158 * t242 + 0.1e1;
t228 = 0.2e1 * (-t242 * t255 + (t187 * t195 * t238 - t211) * t158) / t155 ^ 2;
t227 = 0.2e1 * t255;
t226 = -0.2e1 * t250;
t225 = t171 * t254;
t224 = t158 * t239;
t218 = 0.1e1 + t243;
t217 = t187 * t228;
t216 = -0.2e1 * t187 * t253;
t215 = t253 * t258;
t212 = t174 * t176 * t244;
t210 = t218 * t198;
t209 = t167 * t247 - t249;
t206 = t193 * t233 + t213 * t198;
t170 = -t188 * t231 + t237;
t162 = 0.1e1 / t164;
t161 = t218 * t197 * t176;
t153 = 0.1e1 / t155;
t149 = (t256 * t187 * t173 - t197 * t212) * t198;
t148 = -t188 * t246 + t245 + (t173 * t188 - t174 * t233) * t161;
t146 = -t218 * t215 + (qJD(1) * t210 + t207 * t258) * t176;
t143 = t226 + 0.2e1 * (-t151 * t167 * t162 + (-t162 * t254 - t167 * t250) * t171) * t171;
t142 = t209 * t226 * t239 + (t209 * t188 * t234 + (-t209 * t230 + ((-t166 * t192 - 0.2e1 * t225) * t191 + (-t151 * t191 + (-t171 * t192 + t152) * t190) * t167) * t198) * t187) * t162;
t141 = (t148 * t252 - t157 * t188) * t198 * t228 + ((-t157 * t230 + (-t148 * t193 - t145) * t251) * t188 + (-t157 * t234 - (-t146 * t174 * t197 - t257 * t173 + (t150 * t246 - t173 * t193 - t174 * t229) * t161) * t224 + (t158 * t230 + t198 * t227) * t148 - ((t146 - t229) * t173 + ((-t161 * t197 + 0.1e1) * t193 + (t161 - t197) * t150) * t174) * t188 * t251) * t187) * t153;
t1 = [t198 * t184 * t216 + (t193 * t210 - t230 * t240) * t176, 0, t146, t146, 0, 0; (t157 * t217 + (-t157 * t238 + (qJD(1) * t149 + t145) * t252) * t153) * t197 + (t158 * t217 * t149 + (-((t216 - t238 + (t150 * t184 * t241 + t238) * t176) * t173 + (t215 * t244 - t150 * t187 + (-t182 * t222 + (t150 - 0.2e1 * t235) * t187) * t176) * t174) * t224 + (-t158 * t238 + t187 * t227) * t149 + (-t157 + ((-t194 + t195) * t212 + t256 * t223) * t158) * t187 * qJD(1)) * t153) * t198, 0, t141, t141, 0, 0; 0.2e1 * (t166 * t208 + t170 * t248) * t250 + (0.2e1 * t170 * t225 - t214 * t166 * t231 + t206 * t249 + (-t214 * t171 * t232 + t170 * t151 + t152 * t208 - t206 * t247) * t167) * t162, 0, t142, t142, t143, t143;];
JaD_rot  = t1;
