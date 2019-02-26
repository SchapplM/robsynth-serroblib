% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:26
% EndTime: 2019-02-26 20:16:26
% DurationCPUTime: 0.89s
% Computational Cost: add. (3002->109), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->102)
t207 = sin(pkin(11));
t209 = cos(pkin(11));
t213 = sin(qJ(2));
t210 = cos(pkin(6));
t216 = cos(qJ(2));
t242 = t210 * t216;
t197 = -t207 * t213 + t209 * t242;
t190 = t197 * qJD(2);
t243 = t210 * t213;
t198 = t207 * t216 + t209 * t243;
t212 = sin(qJ(3));
t208 = sin(pkin(6));
t246 = t208 * t212;
t232 = t209 * t246;
t215 = cos(qJ(3));
t239 = qJD(3) * t215;
t162 = -qJD(3) * t232 + t190 * t212 + t198 * t239;
t245 = t208 * t215;
t182 = t198 * t212 + t209 * t245;
t180 = t182 ^ 2;
t201 = -t210 * t215 + t213 * t246;
t195 = 0.1e1 / t201 ^ 2;
t176 = t180 * t195 + 0.1e1;
t174 = 0.1e1 / t176;
t202 = t210 * t212 + t213 * t245;
t240 = qJD(2) * t216;
t231 = t208 * t240;
t187 = t202 * qJD(3) + t212 * t231;
t194 = 0.1e1 / t201;
t250 = t182 * t195;
t146 = (-t162 * t194 + t187 * t250) * t174;
t177 = atan2(-t182, t201);
t172 = sin(t177);
t173 = cos(t177);
t228 = -t172 * t201 - t173 * t182;
t142 = t228 * t146 - t172 * t162 + t173 * t187;
t156 = -t172 * t182 + t173 * t201;
t153 = 0.1e1 / t156;
t154 = 0.1e1 / t156 ^ 2;
t262 = t142 * t153 * t154;
t233 = t207 * t243;
t200 = t209 * t216 - t233;
t225 = -t200 * t212 + t207 * t245;
t261 = -0.2e1 * t225 * t262;
t244 = t208 * t216;
t224 = -t194 * t197 + t244 * t250;
t260 = t212 * t224;
t249 = t187 * t194 * t195;
t259 = -0.2e1 * (t162 * t250 - t180 * t249) / t176 ^ 2;
t186 = t200 * t215 + t207 * t246;
t199 = t207 * t242 + t209 * t213;
t211 = sin(qJ(4));
t214 = cos(qJ(4));
t171 = t186 * t214 + t199 * t211;
t167 = 0.1e1 / t171;
t168 = 0.1e1 / t171 ^ 2;
t258 = t154 * t225;
t192 = t199 * qJD(2);
t165 = t225 * qJD(3) - t192 * t215;
t193 = -qJD(2) * t233 + t209 * t240;
t170 = t186 * t211 - t199 * t214;
t238 = qJD(4) * t170;
t158 = t165 * t214 + t193 * t211 - t238;
t257 = t158 * t167 * t168;
t157 = t171 * qJD(4) + t165 * t211 - t193 * t214;
t166 = t170 ^ 2;
t161 = t166 * t168 + 0.1e1;
t254 = t168 * t170;
t256 = 0.1e1 / t161 ^ 2 * (t157 * t254 - t166 * t257);
t255 = t167 * t211;
t253 = t170 * t214;
t252 = t172 * t225;
t251 = t173 * t225;
t248 = t199 * t212;
t247 = t199 * t215;
t241 = qJD(2) * t213;
t181 = t225 ^ 2;
t152 = t154 * t181 + 0.1e1;
t164 = t186 * qJD(3) - t192 * t212;
t237 = 0.2e1 * (-t164 * t258 - t181 * t262) / t152 ^ 2;
t235 = -0.2e1 * t256;
t234 = t170 * t257;
t230 = -0.2e1 * t182 * t249;
t229 = qJD(4) * t247 - t192;
t227 = t168 * t253 - t255;
t184 = t198 * t215 - t232;
t226 = -t184 * t194 + t202 * t250;
t223 = qJD(3) * t248 + qJD(4) * t200 - t193 * t215;
t191 = t198 * qJD(2);
t188 = -t201 * qJD(3) + t215 * t231;
t179 = t200 * t211 - t214 * t247;
t178 = -t200 * t214 - t211 * t247;
t163 = -t182 * qJD(3) + t190 * t215;
t159 = 0.1e1 / t161;
t149 = 0.1e1 / t152;
t148 = t174 * t260;
t147 = t226 * t174;
t144 = (-t172 * t197 + t173 * t244) * t212 + t228 * t148;
t143 = t228 * t147 - t172 * t184 + t173 * t202;
t141 = t226 * t259 + (t202 * t230 - t163 * t194 + (t162 * t202 + t182 * t188 + t184 * t187) * t195) * t174;
t139 = t259 * t260 + (t224 * t239 + (t230 * t244 + t191 * t194 + (t187 * t197 + (t162 * t216 - t182 * t241) * t208) * t195) * t212) * t174;
t1 = [0, t139, t141, 0, 0, 0; 0 (-t144 * t258 + t153 * t248) * t237 + ((-t193 * t212 - t199 * t239) * t153 + t144 * t261 + (-t144 * t164 + t248 * t142 + (-t139 * t182 - t148 * t162 + (-t212 * t241 + t216 * t239) * t208 + (-t148 * t201 - t197 * t212) * t146) * t251 + (-t197 * t239 - t139 * t201 - t148 * t187 + t191 * t212 + (t148 * t182 - t212 * t244) * t146) * t252) * t154) * t149 (-t143 * t258 - t153 * t186) * t237 + (t143 * t261 + t165 * t153 + (-t186 * t142 - t143 * t164 + (-t141 * t182 - t147 * t162 + t188 + (-t147 * t201 - t184) * t146) * t251 + (-t141 * t201 - t147 * t187 - t163 + (t147 * t182 - t202) * t146) * t252) * t154) * t149, 0, 0, 0; 0, 0.2e1 * (-t167 * t178 + t179 * t254) * t256 + (0.2e1 * t179 * t234 - t229 * t167 * t214 + t223 * t255 + (-t229 * t170 * t211 - t179 * t157 - t178 * t158 - t223 * t253) * t168) * t159, -t227 * t225 * t235 + (t227 * t164 - ((-qJD(4) * t167 - 0.2e1 * t234) * t214 + (t157 * t214 + (t158 - t238) * t211) * t168) * t225) * t159, t235 + 0.2e1 * (t157 * t168 * t159 + (-t159 * t257 - t168 * t256) * t170) * t170, 0, 0;];
JaD_rot  = t1;
