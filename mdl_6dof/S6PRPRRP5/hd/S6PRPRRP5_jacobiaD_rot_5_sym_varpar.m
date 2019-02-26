% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:31
% EndTime: 2019-02-26 19:52:32
% DurationCPUTime: 0.90s
% Computational Cost: add. (3002->108), mult. (9085->231), div. (559->12), fcn. (11668->13), ass. (0->104)
t208 = sin(pkin(10));
t210 = cos(pkin(10));
t217 = cos(qJ(2));
t211 = cos(pkin(6));
t214 = sin(qJ(2));
t243 = t211 * t214;
t202 = t208 * t217 + t210 * t243;
t196 = t202 * qJD(2);
t209 = sin(pkin(6));
t216 = cos(qJ(4));
t213 = sin(qJ(4));
t242 = t211 * t217;
t227 = -t208 * t214 + t210 * t242;
t225 = t227 * t213;
t167 = qJD(4) * t225 + (qJD(4) * t209 * t210 + t196) * t216;
t247 = t209 * t213;
t186 = t210 * t247 - t227 * t216;
t183 = t186 ^ 2;
t244 = t209 * t217;
t205 = t211 * t213 + t216 * t244;
t200 = 0.1e1 / t205 ^ 2;
t178 = t183 * t200 + 0.1e1;
t176 = 0.1e1 / t178;
t206 = t211 * t216 - t213 * t244;
t246 = t209 * t214;
t232 = qJD(2) * t246;
t189 = t206 * qJD(4) - t216 * t232;
t199 = 0.1e1 / t205;
t253 = t186 * t200;
t148 = (t167 * t199 - t189 * t253) * t176;
t179 = atan2(t186, t205);
t174 = sin(t179);
t175 = cos(t179);
t230 = -t174 * t205 + t175 * t186;
t144 = t230 * t148 + t174 * t167 + t175 * t189;
t158 = t174 * t186 + t175 * t205;
t155 = 0.1e1 / t158;
t156 = 0.1e1 / t158 ^ 2;
t265 = t144 * t155 * t156;
t203 = t208 * t242 + t210 * t214;
t184 = -t203 * t216 + t208 * t247;
t264 = 0.2e1 * t184 * t265;
t251 = t189 * t199 * t200;
t263 = (t167 * t253 - t183 * t251) / t178 ^ 2;
t234 = t186 * t246;
t226 = t199 * t202 + t200 * t234;
t262 = t216 * t226;
t245 = t209 * t216;
t185 = t203 * t213 + t208 * t245;
t233 = t208 * t243;
t204 = t210 * t217 - t233;
t212 = sin(qJ(5));
t215 = cos(qJ(5));
t173 = t185 * t215 + t204 * t212;
t169 = 0.1e1 / t173;
t170 = 0.1e1 / t173 ^ 2;
t261 = t156 * t184;
t241 = qJD(2) * t217;
t198 = -qJD(2) * t233 + t210 * t241;
t164 = -t184 * qJD(4) + t198 * t213;
t197 = t203 * qJD(2);
t249 = t204 * t215;
t172 = t185 * t212 - t249;
t239 = qJD(5) * t172;
t160 = t164 * t215 - t197 * t212 - t239;
t260 = t160 * t169 * t170;
t159 = t173 * qJD(5) + t164 * t212 + t197 * t215;
t168 = t172 ^ 2;
t163 = t168 * t170 + 0.1e1;
t257 = t170 * t172;
t259 = 0.1e1 / t163 ^ 2 * (t159 * t257 - t168 * t260);
t258 = t169 * t212;
t256 = t172 * t215;
t255 = t174 * t184;
t254 = t175 * t184;
t252 = t186 * t206;
t250 = t204 * t213;
t248 = t204 * t216;
t240 = qJD(4) * t213;
t182 = t184 ^ 2;
t154 = t156 * t182 + 0.1e1;
t165 = t185 * qJD(4) - t198 * t216;
t238 = 0.2e1 * (t165 * t261 - t182 * t265) / t154 ^ 2;
t236 = -0.2e1 * t259;
t235 = t172 * t260;
t231 = qJD(5) * t250 + t198;
t229 = t170 * t256 - t258;
t187 = t210 * t245 + t225;
t228 = -t187 * t199 + t200 * t252;
t224 = qJD(4) * t248 - qJD(5) * t203 - t197 * t213;
t195 = t227 * qJD(2);
t188 = -t205 * qJD(4) + t213 * t232;
t181 = -t203 * t212 + t213 * t249;
t180 = t203 * t215 + t212 * t250;
t166 = t186 * qJD(4) + t196 * t213;
t161 = 0.1e1 / t163;
t151 = 0.1e1 / t154;
t150 = t176 * t262;
t149 = t228 * t176;
t146 = (t174 * t202 - t175 * t246) * t216 + t230 * t150;
t145 = -t230 * t149 + t174 * t187 + t175 * t206;
t143 = 0.2e1 * t228 * t263 + (0.2e1 * t251 * t252 - t166 * t199 + (-t167 * t206 - t186 * t188 - t187 * t189) * t200) * t176;
t141 = -0.2e1 * t262 * t263 + (-t226 * t240 + (-0.2e1 * t234 * t251 + t195 * t199 + (-t189 * t202 + (t167 * t214 + t186 * t241) * t209) * t200) * t216) * t176;
t1 = [0, t141, 0, t143, 0, 0; 0 (t146 * t261 + t155 * t248) * t238 + ((t197 * t216 + t204 * t240) * t155 + t146 * t264 + (-t146 * t165 + t248 * t144 - (t141 * t186 + t150 * t167 + (t214 * t240 - t216 * t241) * t209 + (-t150 * t205 + t202 * t216) * t148) * t254 - (-t202 * t240 - t141 * t205 - t150 * t189 + t195 * t216 + (-t150 * t186 + t214 * t245) * t148) * t255) * t156) * t151, 0 (t145 * t261 - t155 * t185) * t238 + (t145 * t264 + t164 * t155 + (-t185 * t144 - t145 * t165 - (t143 * t186 - t149 * t167 + t188 + (t149 * t205 + t187) * t148) * t254 - (-t143 * t205 + t149 * t189 - t166 + (t149 * t186 - t206) * t148) * t255) * t156) * t151, 0, 0; 0, 0.2e1 * (-t169 * t180 + t181 * t257) * t259 + (0.2e1 * t181 * t235 + t231 * t169 * t215 + t224 * t258 + (t231 * t172 * t212 - t181 * t159 - t180 * t160 - t224 * t256) * t170) * t161, 0, t229 * t184 * t236 + (t229 * t165 + ((-qJD(5) * t169 - 0.2e1 * t235) * t215 + (t159 * t215 + (t160 - t239) * t212) * t170) * t184) * t161, t236 + 0.2e1 * (t159 * t170 * t161 + (-t161 * t260 - t170 * t259) * t172) * t172, 0;];
JaD_rot  = t1;
