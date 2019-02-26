% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:17
% EndTime: 2019-02-26 22:18:18
% DurationCPUTime: 0.76s
% Computational Cost: add. (4138->97), mult. (4025->207), div. (771->12), fcn. (4686->9), ass. (0->100)
t208 = sin(qJ(1));
t207 = qJ(2) + qJ(3);
t199 = sin(t207);
t194 = 0.1e1 / t199 ^ 2;
t201 = cos(t207);
t197 = t201 ^ 2;
t254 = t194 * t197;
t229 = 0.1e1 + t254;
t271 = t208 * t229;
t202 = qJD(5) + qJD(6);
t225 = qJD(1) * t199 + t202;
t203 = qJD(2) + qJD(3);
t209 = cos(qJ(1));
t246 = t203 * t209;
t270 = -t201 * t246 + t208 * t225;
t244 = t208 * t201;
t269 = t203 * t244 + t209 * t225;
t188 = atan2(-t244, t199);
t187 = cos(t188);
t186 = sin(t188);
t235 = t186 * t244;
t173 = t187 * t199 - t235;
t170 = 0.1e1 / t173;
t206 = qJ(5) + qJ(6);
t198 = sin(t206);
t200 = cos(t206);
t245 = t208 * t200;
t249 = t199 * t209;
t183 = t198 * t249 + t245;
t179 = 0.1e1 / t183;
t193 = 0.1e1 / t199;
t171 = 0.1e1 / t173 ^ 2;
t180 = 0.1e1 / t183 ^ 2;
t204 = t208 ^ 2;
t191 = t204 * t254 + 0.1e1;
t189 = 0.1e1 / t191;
t268 = t189 - 0.1e1;
t242 = qJD(1) * t209;
t230 = t201 * t242;
t247 = t203 * t208;
t233 = t194 * t247;
t163 = ((t247 * t199 - t230) * t193 + t197 * t233) * t189;
t256 = t187 * t201;
t158 = (-t163 * t208 + t203) * t256 + (-t230 + (-t163 + t247) * t199) * t186;
t267 = t158 * t170 * t171;
t266 = t163 * t186;
t226 = t199 * t202 + qJD(1);
t220 = t226 * t209;
t165 = -t270 * t198 + t200 * t220;
t265 = t165 * t179 * t180;
t196 = t201 * t197;
t255 = t193 * t201;
t218 = t203 * (-t193 * t194 * t196 - t255);
t252 = t197 * t208;
t223 = t242 * t252;
t264 = (t194 * t223 + t204 * t218) / t191 ^ 2;
t263 = t171 * t201;
t262 = t171 * t209;
t164 = t198 * t220 + t270 * t200;
t248 = t200 * t209;
t182 = t198 * t208 - t199 * t248;
t178 = t182 ^ 2;
t177 = t178 * t180 + 0.1e1;
t258 = t180 * t182;
t261 = 0.1e1 / t177 ^ 2 * (t164 * t258 - t178 * t265);
t176 = t189 * t271;
t260 = t176 * t208;
t259 = t179 * t200;
t257 = t182 * t198;
t205 = t209 ^ 2;
t253 = t197 * t205;
t251 = t199 * t203;
t250 = t199 * t208;
t243 = qJD(1) * t208;
t168 = t171 * t253 + 0.1e1;
t241 = 0.2e1 * (-t253 * t267 + (-t201 * t205 * t251 - t223) * t171) / t168 ^ 2;
t240 = 0.2e1 * t267;
t239 = -0.2e1 * t264;
t238 = 0.2e1 * t261;
t237 = t201 * t264;
t236 = t201 * t262;
t234 = t193 * t252;
t228 = t201 * t241;
t227 = 0.2e1 * t182 * t265;
t224 = t187 * t189 * t193 * t197;
t222 = t229 * t209;
t221 = t226 * t208;
t219 = t257 * t180 + t259;
t217 = t219 * t209;
t185 = -t198 * t250 + t248;
t184 = t198 * t209 + t199 * t245;
t174 = 0.1e1 / t177;
t166 = 0.1e1 / t168;
t162 = (t186 * t201 * t268 + t208 * t224) * t209;
t161 = t186 * t250 + t256 + (-t186 * t199 - t187 * t244) * t176;
t159 = t239 * t271 + (qJD(1) * t222 + 0.2e1 * t208 * t218) * t189;
t156 = -0.2e1 * t261 + 0.2e1 * (t164 * t180 * t174 + (-t174 * t265 - t180 * t261) * t182) * t182;
t155 = t201 * t217 * t238 + (t217 * t251 + (t219 * t243 + ((t179 * t202 + t227) * t198 + (-t164 * t198 + (-t182 * t202 + t165) * t200) * t180) * t209) * t201) * t174;
t154 = (t161 * t263 + t170 * t199) * t209 * t241 + ((t170 * t243 + (t161 * t203 + t158) * t262) * t199 + (-t170 * t246 - (-t159 * t187 * t208 + t186 * t247 + t260 * t266 - t266 + (-t186 * t203 - t187 * t242) * t176) * t236 + (t171 * t243 + t209 * t240) * t161 - ((-t159 + t242) * t186 + ((-0.1e1 + t260) * t203 + (-t176 + t208) * t163) * t187) * t171 * t249) * t201) * t166;
t1 = [0.2e1 * t209 * t193 * t237 + (t203 * t222 + t243 * t255) * t189, t159, t159, 0, 0, 0; (t170 * t228 + (t170 * t251 + (qJD(1) * t162 + t158) * t263) * t166) * t208 + (t171 * t228 * t162 + (-((-0.2e1 * t237 + t251 + (-t163 * t234 - t251) * t189) * t186 + (t234 * t239 - t163 * t201 + (-t196 * t233 + (t163 - 0.2e1 * t247) * t201) * t189) * t187) * t236 + (t171 * t251 + t201 * t240) * t162 + (-t170 + ((t204 - t205) * t224 + t268 * t235) * t171) * t201 * qJD(1)) * t166) * t209, t154, t154, 0, 0, 0; (-t179 * t184 + t185 * t258) * t238 + (t185 * t227 - t179 * t198 * t221 + t269 * t259 + (t182 * t200 * t221 - t185 * t164 - t184 * t165 + t269 * t257) * t180) * t174, t155, t155, 0, t156, t156;];
JaD_rot  = t1;
