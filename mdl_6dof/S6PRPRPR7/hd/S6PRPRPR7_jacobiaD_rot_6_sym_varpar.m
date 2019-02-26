% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:47
% EndTime: 2019-02-26 19:49:48
% DurationCPUTime: 0.93s
% Computational Cost: add. (3002->108), mult. (9085->223), div. (559->12), fcn. (11668->13), ass. (0->102)
t207 = cos(pkin(10));
t213 = cos(qJ(4));
t205 = sin(pkin(10));
t211 = sin(qJ(2));
t208 = cos(pkin(6));
t214 = cos(qJ(2));
t240 = t208 * t214;
t224 = -t205 * t211 + t207 * t240;
t206 = sin(pkin(6));
t210 = sin(qJ(4));
t245 = t206 * t210;
t186 = t207 * t245 - t224 * t213;
t241 = t208 * t211;
t199 = t205 * t214 + t207 * t241;
t193 = t199 * qJD(2);
t165 = t186 * qJD(4) + t193 * t210;
t243 = t206 * t213;
t187 = t207 * t243 + t224 * t210;
t183 = t187 ^ 2;
t242 = t206 * t214;
t223 = -t208 * t213 + t210 * t242;
t197 = 0.1e1 / t223 ^ 2;
t177 = t183 * t197 + 0.1e1;
t175 = 0.1e1 / t177;
t202 = -t208 * t210 - t213 * t242;
t244 = t206 * t211;
t232 = qJD(2) * t244;
t188 = t202 * qJD(4) + t210 * t232;
t196 = 0.1e1 / t223;
t249 = t187 * t197;
t147 = (t165 * t196 - t188 * t249) * t175;
t178 = atan2(t187, -t223);
t173 = sin(t178);
t174 = cos(t178);
t228 = t173 * t223 + t174 * t187;
t143 = t228 * t147 - t173 * t165 + t174 * t188;
t157 = t173 * t187 - t174 * t223;
t154 = 0.1e1 / t157;
t155 = 0.1e1 / t157 ^ 2;
t262 = t143 * t154 * t155;
t200 = t205 * t240 + t207 * t211;
t184 = -t200 * t213 + t205 * t245;
t233 = t205 * t241;
t201 = t207 * t214 - t233;
t209 = sin(qJ(6));
t212 = cos(qJ(6));
t227 = t184 * t212 - t201 * t209;
t261 = t227 * qJD(6);
t185 = t200 * t210 + t205 * t243;
t260 = 0.2e1 * t185 * t262;
t222 = -t196 * t199 + t244 * t249;
t259 = t210 * t222;
t172 = t184 * t209 + t201 * t212;
t168 = 0.1e1 / t172;
t169 = 0.1e1 / t172 ^ 2;
t239 = qJD(2) * t214;
t195 = -qJD(2) * t233 + t207 * t239;
t164 = t185 * qJD(4) - t195 * t213;
t194 = t200 * qJD(2);
t158 = t172 * qJD(6) - t164 * t212 - t194 * t209;
t167 = t227 ^ 2;
t162 = t167 * t169 + 0.1e1;
t253 = t169 * t227;
t159 = t164 * t209 - t194 * t212 + t261;
t256 = t159 * t168 * t169;
t258 = (-t158 * t253 - t167 * t256) / t162 ^ 2;
t257 = t155 * t185;
t163 = -t184 * qJD(4) + t195 * t210;
t255 = t163 * t155;
t254 = t168 * t212;
t252 = t227 * t209;
t251 = t173 * t185;
t250 = t174 * t185;
t248 = t188 * t196 * t197;
t247 = t201 * t210;
t246 = t201 * t213;
t238 = qJD(4) * t213;
t182 = t185 ^ 2;
t153 = t182 * t155 + 0.1e1;
t237 = 0.2e1 * (-t182 * t262 + t185 * t255) / t153 ^ 2;
t236 = 0.2e1 * t258;
t235 = 0.2e1 * (-t165 * t249 + t183 * t248) / t177 ^ 2;
t231 = -0.2e1 * t227 * t256;
t230 = -0.2e1 * t187 * t248;
t229 = -qJD(6) * t246 - t195;
t226 = -t169 * t252 + t254;
t225 = -t186 * t196 + t202 * t249;
t221 = qJD(4) * t247 + qJD(6) * t200 + t194 * t213;
t192 = t224 * qJD(2);
t189 = t223 * qJD(4) + t213 * t232;
t180 = -t200 * t212 - t209 * t246;
t179 = -t200 * t209 + t212 * t246;
t166 = t187 * qJD(4) + t193 * t213;
t160 = 0.1e1 / t162;
t150 = 0.1e1 / t153;
t149 = t175 * t259;
t148 = t225 * t175;
t145 = (-t173 * t199 + t174 * t244) * t210 - t228 * t149;
t144 = -t228 * t148 - t173 * t186 + t174 * t202;
t142 = t225 * t235 + (t202 * t230 + t166 * t196 + (t165 * t202 + t186 * t188 - t187 * t189) * t197) * t175;
t140 = t235 * t259 + (-t222 * t238 + (t230 * t244 + t192 * t196 + (t188 * t199 + (t165 * t211 - t187 * t239) * t206) * t197) * t210) * t175;
t1 = [0, t140, 0, t142, 0, 0; 0 (t145 * t257 - t154 * t247) * t237 + ((-t194 * t210 + t201 * t238) * t154 + (-t255 + t260) * t145 + (-t247 * t143 - (t140 * t187 + t149 * t165 + (t210 * t239 + t211 * t238) * t206 + (-t149 * t223 - t199 * t210) * t147) * t250 - (-t199 * t238 + t140 * t223 + t149 * t188 - t192 * t210 + (t149 * t187 - t210 * t244) * t147) * t251) * t155) * t150, 0 (t144 * t257 + t154 * t184) * t237 + (t144 * t260 - t164 * t154 + (t184 * t143 - t144 * t163 - (t142 * t187 + t148 * t165 + t189 + (-t148 * t223 - t186) * t147) * t250 - (t142 * t223 + t148 * t188 - t166 + (t148 * t187 - t202) * t147) * t251) * t155) * t150, 0, 0; 0 (-t168 * t179 - t180 * t253) * t236 + (t180 * t231 + t229 * t168 * t209 - t221 * t254 + (t212 * t227 * t229 - t180 * t158 - t179 * t159 + t221 * t252) * t169) * t160, 0, t226 * t185 * t236 + (-t226 * t163 + ((qJD(6) * t168 + t231) * t209 + (-t158 * t209 + (t159 + t261) * t212) * t169) * t185) * t160, 0, -0.2e1 * t258 - 0.2e1 * (t158 * t169 * t160 - (-t160 * t256 - t169 * t258) * t227) * t227;];
JaD_rot  = t1;
