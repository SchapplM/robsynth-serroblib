% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_rot = S6PRPRRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:31
% EndTime: 2019-02-26 19:52:32
% DurationCPUTime: 0.88s
% Computational Cost: add. (3002->108), mult. (9085->231), div. (559->12), fcn. (11668->13), ass. (0->105)
t218 = sin(pkin(10));
t220 = cos(pkin(10));
t227 = cos(qJ(2));
t221 = cos(pkin(6));
t224 = sin(qJ(2));
t253 = t221 * t224;
t212 = t218 * t227 + t220 * t253;
t206 = t212 * qJD(2);
t219 = sin(pkin(6));
t226 = cos(qJ(4));
t223 = sin(qJ(4));
t252 = t221 * t227;
t237 = -t218 * t224 + t220 * t252;
t235 = t237 * t223;
t177 = qJD(4) * t235 + (qJD(4) * t219 * t220 + t206) * t226;
t257 = t219 * t223;
t196 = t220 * t257 - t237 * t226;
t193 = t196 ^ 2;
t254 = t219 * t227;
t215 = t221 * t223 + t226 * t254;
t210 = 0.1e1 / t215 ^ 2;
t188 = t193 * t210 + 0.1e1;
t186 = 0.1e1 / t188;
t216 = t221 * t226 - t223 * t254;
t256 = t219 * t224;
t242 = qJD(2) * t256;
t199 = t216 * qJD(4) - t226 * t242;
t209 = 0.1e1 / t215;
t263 = t196 * t210;
t158 = (t177 * t209 - t199 * t263) * t186;
t189 = atan2(t196, t215);
t184 = sin(t189);
t185 = cos(t189);
t240 = -t184 * t215 + t185 * t196;
t154 = t240 * t158 + t184 * t177 + t185 * t199;
t168 = t184 * t196 + t185 * t215;
t165 = 0.1e1 / t168;
t166 = 0.1e1 / t168 ^ 2;
t276 = t154 * t165 * t166;
t213 = t218 * t252 + t220 * t224;
t194 = -t213 * t226 + t218 * t257;
t275 = 0.2e1 * t194 * t276;
t261 = t199 * t209 * t210;
t274 = (t177 * t263 - t193 * t261) / t188 ^ 2;
t244 = t196 * t256;
t236 = t209 * t212 + t210 * t244;
t273 = t226 * t236;
t255 = t219 * t226;
t195 = t213 * t223 + t218 * t255;
t243 = t218 * t253;
t214 = t220 * t227 - t243;
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t183 = t195 * t225 + t214 * t222;
t179 = 0.1e1 / t183;
t180 = 0.1e1 / t183 ^ 2;
t251 = qJD(2) * t227;
t208 = -qJD(2) * t243 + t220 * t251;
t174 = -t194 * qJD(4) + t208 * t223;
t207 = t213 * qJD(2);
t169 = t183 * qJD(5) + t174 * t222 + t207 * t225;
t259 = t214 * t225;
t182 = t195 * t222 - t259;
t178 = t182 ^ 2;
t173 = t178 * t180 + 0.1e1;
t267 = t180 * t182;
t249 = qJD(5) * t182;
t170 = t174 * t225 - t207 * t222 - t249;
t270 = t170 * t179 * t180;
t272 = (t169 * t267 - t178 * t270) / t173 ^ 2;
t271 = t166 * t194;
t175 = t195 * qJD(4) - t208 * t226;
t269 = t175 * t166;
t268 = t179 * t222;
t266 = t182 * t225;
t265 = t184 * t194;
t264 = t185 * t194;
t262 = t196 * t216;
t260 = t214 * t223;
t258 = t214 * t226;
t250 = qJD(4) * t223;
t192 = t194 ^ 2;
t164 = t166 * t192 + 0.1e1;
t248 = 0.2e1 * (-t192 * t276 + t194 * t269) / t164 ^ 2;
t247 = -0.2e1 * t272;
t245 = t182 * t270;
t241 = qJD(5) * t260 + t208;
t239 = t180 * t266 - t268;
t197 = t220 * t255 + t235;
t238 = -t197 * t209 + t210 * t262;
t234 = qJD(4) * t258 - qJD(5) * t213 - t207 * t223;
t205 = t237 * qJD(2);
t198 = -t215 * qJD(4) + t223 * t242;
t191 = -t213 * t222 + t223 * t259;
t190 = t213 * t225 + t222 * t260;
t176 = t196 * qJD(4) + t206 * t223;
t171 = 0.1e1 / t173;
t161 = 0.1e1 / t164;
t160 = t186 * t273;
t159 = t238 * t186;
t156 = (t184 * t212 - t185 * t256) * t226 + t240 * t160;
t155 = -t240 * t159 + t184 * t197 + t185 * t216;
t153 = 0.2e1 * t238 * t274 + (0.2e1 * t261 * t262 - t176 * t209 + (-t177 * t216 - t196 * t198 - t197 * t199) * t210) * t186;
t151 = -0.2e1 * t273 * t274 + (-t236 * t250 + (-0.2e1 * t244 * t261 + t205 * t209 + (-t199 * t212 + (t177 * t224 + t196 * t251) * t219) * t210) * t226) * t186;
t1 = [0, t151, 0, t153, 0, 0; 0 (t156 * t271 + t165 * t258) * t248 + ((t207 * t226 + t214 * t250) * t165 + (-t269 + t275) * t156 + (t258 * t154 - (t151 * t196 + t160 * t177 + (t224 * t250 - t226 * t251) * t219 + (-t160 * t215 + t212 * t226) * t158) * t264 - (-t212 * t250 - t151 * t215 - t160 * t199 + t205 * t226 + (-t160 * t196 + t224 * t255) * t158) * t265) * t166) * t161, 0 (t155 * t271 - t165 * t195) * t248 + (t155 * t275 + t174 * t165 + (-t195 * t154 - t155 * t175 - (t153 * t196 - t159 * t177 + t198 + (t159 * t215 + t197) * t158) * t264 - (-t153 * t215 + t159 * t199 - t176 + (t159 * t196 - t216) * t158) * t265) * t166) * t161, 0, 0; 0, 0.2e1 * (-t179 * t190 + t191 * t267) * t272 + (0.2e1 * t191 * t245 + t241 * t179 * t225 + t234 * t268 + (t241 * t182 * t222 - t191 * t169 - t190 * t170 - t234 * t266) * t180) * t171, 0, t239 * t194 * t247 + (t239 * t175 + ((-qJD(5) * t179 - 0.2e1 * t245) * t225 + (t169 * t225 + (t170 - t249) * t222) * t180) * t194) * t171, t247 + 0.2e1 * (t169 * t171 * t180 + (-t171 * t270 - t180 * t272) * t182) * t182, 0;];
JaD_rot  = t1;
