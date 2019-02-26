% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:13
% EndTime: 2019-02-26 21:41:14
% DurationCPUTime: 1.20s
% Computational Cost: add. (7491->94), mult. (9702->187), div. (711->12), fcn. (12089->11), ass. (0->95)
t274 = qJ(4) + pkin(10);
t263 = sin(t274);
t264 = cos(t274);
t293 = sin(qJ(2));
t295 = cos(qJ(2));
t211 = -t263 * t295 + t264 * t293;
t294 = sin(qJ(1));
t202 = t211 * t294;
t210 = t263 * t293 + t264 * t295;
t189 = atan2(t202, t210);
t184 = sin(t189);
t185 = cos(t189);
t204 = t210 * t294;
t258 = -t184 * t210 + t185 * t202;
t200 = t202 ^ 2;
t208 = 0.1e1 / t210 ^ 2;
t188 = t200 * t208 + 0.1e1;
t186 = 0.1e1 / t188;
t207 = 0.1e1 / t210;
t276 = t202 * t211;
t250 = -t204 * t207 - t208 * t276;
t303 = t250 * t186;
t161 = -t184 * t204 + t185 * t211 + t258 * t303;
t296 = cos(qJ(1));
t206 = t211 * t296;
t322 = t161 * t206;
t174 = t184 * t202 + t185 * t210;
t171 = 0.1e1 / t174;
t172 = 0.1e1 / t174 ^ 2;
t205 = t210 * t296;
t201 = t206 ^ 2;
t170 = t172 * t201 + 0.1e1;
t311 = qJD(2) - qJD(4);
t178 = -t202 * qJD(1) + t205 * t311;
t284 = t178 * t172;
t179 = -t206 * qJD(1) - t204 * t311;
t191 = t311 * t211;
t277 = t202 * t208;
t252 = -t179 * t207 + t191 * t277;
t164 = t252 * t186;
t159 = t164 * t258 - t184 * t179 - t185 * t191;
t292 = t159 * t171 * t172;
t273 = 0.2e1 * (-t201 * t292 + t206 * t284) / t170 ^ 2;
t321 = (t171 * t205 + t172 * t322) * t273;
t177 = t204 * qJD(1) + t206 * t311;
t229 = sin(qJ(6));
t230 = cos(qJ(6));
t199 = t205 * t230 - t229 * t294;
t266 = qJD(1) * t296;
t175 = qJD(6) * t199 - t177 * t229 + t230 * t266;
t248 = -t205 * t229 - t230 * t294;
t312 = t248 * qJD(6);
t176 = -t177 * t230 - t229 * t266 + t312;
t192 = t248 ^ 2;
t194 = 0.1e1 / t199 ^ 2;
t183 = t192 * t194 + 0.1e1;
t181 = 0.1e1 / t183;
t193 = 0.1e1 / t199;
t279 = t194 * t248;
t251 = -t193 * t229 - t230 * t279;
t286 = t176 * t193 * t194;
t298 = -0.2e1 * t248;
t262 = t286 * t298;
t320 = (t251 * t178 - (t194 * ((-t176 - t312) * t229 - t175 * t230) + (qJD(6) * t193 + t262) * t230) * t206) * t181;
t319 = -t205 * t159 + t161 * t178;
t272 = 0.2e1 * t292;
t318 = -t177 * t171 - t272 * t322;
t190 = t311 * t210;
t316 = (t210 * t303 + t204) * t164 + t303 * t179 - t190;
t180 = t205 * qJD(1) - t202 * t311;
t315 = -(-t202 * t303 - t211) * t164 - t303 * t191 + t180;
t305 = t191 * t208;
t280 = t207 * t305;
t313 = t186 * ((t179 * t211 - t190 * t202 - t191 * t204) * t208 - t180 * t207 - 0.2e1 * t276 * t280);
t278 = t202 * t207;
t301 = (-t185 * t278 + t184) * t186 - t184;
t297 = -0.2e1 * t206;
t291 = (-t175 * t279 - t192 * t286) / t183 ^ 2;
t290 = (-t179 * t277 + t200 * t280) / t188 ^ 2;
t283 = t181 * t194;
t282 = t184 * t206;
t281 = t185 * t206;
t271 = -0.2e1 * t291;
t270 = -0.2e1 * t290;
t269 = t194 * t291;
t268 = t207 * t290;
t267 = t175 * t283;
t265 = qJD(1) * t294;
t197 = -t204 * t230 - t229 * t296;
t249 = t204 * t229 - t230 * t296;
t168 = 0.1e1 / t170;
t163 = t301 * t206;
t157 = t250 * t270 + t313;
t156 = 0.2e1 * t250 * t290 - t313;
t1 = [t268 * t297 + (t178 * t207 + t206 * t305) * t186, t156, 0, t157, 0, 0; -t202 * t171 * t273 + (-t179 * t171 + (-t159 * t202 - t163 * t178) * t172) * t168 - ((-t163 * t272 + t301 * t284) * t168 + (-t163 * t273 + ((t164 * t186 * t278 + t270) * t282 + (0.2e1 * t202 * t268 - t164 + (t164 - t252) * t186) * t281) * t168) * t172) * t206, t321 + (((t156 * t202 + t316) * t281 + (-t156 * t210 + t315) * t282 - t319) * t172 - t318) * t168, 0, -t321 + (((t157 * t202 - t316) * t281 + (-t157 * t210 - t315) * t282 + t319) * t172 + t318) * t168, 0, 0; (t269 * t298 - t267) * t197 - (-t176 * t283 + t193 * t271) * t249 + ((qJD(6) * t197 - t180 * t229 - t230 * t265) * t193 + (qJD(6) * t249 - t180 * t230 + t229 * t265) * t279 + t197 * t262) * t181, t251 * t291 * t297 + t320, 0, -t251 * t206 * t271 - t320, 0, t271 + (t267 - (-t181 * t286 - t269) * t248) * t298;];
JaD_rot  = t1;
