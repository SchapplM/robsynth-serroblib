% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:25
% EndTime: 2019-02-26 21:31:26
% DurationCPUTime: 1.22s
% Computational Cost: add. (7491->94), mult. (9702->187), div. (711->12), fcn. (12089->11), ass. (0->95)
t272 = pkin(10) + qJ(5);
t261 = sin(t272);
t262 = cos(t272);
t291 = sin(qJ(2));
t293 = cos(qJ(2));
t209 = -t293 * t261 + t291 * t262;
t292 = sin(qJ(1));
t200 = t209 * t292;
t208 = t291 * t261 + t293 * t262;
t187 = atan2(t200, t208);
t182 = sin(t187);
t183 = cos(t187);
t202 = t208 * t292;
t256 = -t182 * t208 + t183 * t200;
t198 = t200 ^ 2;
t206 = 0.1e1 / t208 ^ 2;
t186 = t198 * t206 + 0.1e1;
t184 = 0.1e1 / t186;
t205 = 0.1e1 / t208;
t274 = t200 * t209;
t248 = -t202 * t205 - t206 * t274;
t301 = t248 * t184;
t159 = -t182 * t202 + t183 * t209 + t256 * t301;
t294 = cos(qJ(1));
t204 = t209 * t294;
t320 = t159 * t204;
t172 = t182 * t200 + t183 * t208;
t169 = 0.1e1 / t172;
t170 = 0.1e1 / t172 ^ 2;
t203 = t208 * t294;
t199 = t204 ^ 2;
t168 = t199 * t170 + 0.1e1;
t309 = qJD(2) - qJD(5);
t176 = -t200 * qJD(1) + t309 * t203;
t282 = t176 * t170;
t177 = -t204 * qJD(1) - t309 * t202;
t189 = t309 * t209;
t275 = t200 * t206;
t250 = -t177 * t205 + t189 * t275;
t162 = t250 * t184;
t157 = t162 * t256 - t182 * t177 - t183 * t189;
t290 = t157 * t169 * t170;
t271 = 0.2e1 * (-t199 * t290 + t204 * t282) / t168 ^ 2;
t319 = (t169 * t203 + t170 * t320) * t271;
t175 = qJD(1) * t202 + t309 * t204;
t227 = sin(qJ(6));
t228 = cos(qJ(6));
t197 = t203 * t228 - t227 * t292;
t264 = qJD(1) * t294;
t173 = qJD(6) * t197 - t175 * t227 + t228 * t264;
t246 = -t203 * t227 - t228 * t292;
t310 = qJD(6) * t246;
t174 = -t175 * t228 - t227 * t264 + t310;
t190 = t246 ^ 2;
t192 = 0.1e1 / t197 ^ 2;
t181 = t190 * t192 + 0.1e1;
t179 = 0.1e1 / t181;
t191 = 0.1e1 / t197;
t277 = t192 * t246;
t249 = -t227 * t191 - t228 * t277;
t284 = t174 * t191 * t192;
t296 = -0.2e1 * t246;
t260 = t284 * t296;
t318 = (t249 * t176 - ((t227 * (-t174 - t310) - t173 * t228) * t192 + (qJD(6) * t191 + t260) * t228) * t204) * t179;
t317 = -t203 * t157 + t159 * t176;
t270 = 0.2e1 * t290;
t316 = -t175 * t169 - t270 * t320;
t188 = t309 * t208;
t314 = (t208 * t301 + t202) * t162 + t301 * t177 - t188;
t178 = qJD(1) * t203 - t309 * t200;
t313 = -(-t200 * t301 - t209) * t162 - t301 * t189 + t178;
t303 = t189 * t206;
t278 = t205 * t303;
t311 = ((t177 * t209 - t188 * t200 - t189 * t202) * t206 - t178 * t205 - 0.2e1 * t274 * t278) * t184;
t276 = t200 * t205;
t299 = t184 * (-t183 * t276 + t182) - t182;
t295 = -0.2e1 * t204;
t289 = (-t173 * t277 - t190 * t284) / t181 ^ 2;
t288 = (-t177 * t275 + t198 * t278) / t186 ^ 2;
t281 = t179 * t192;
t280 = t182 * t204;
t279 = t183 * t204;
t269 = -0.2e1 * t289;
t268 = -0.2e1 * t288;
t267 = t192 * t289;
t266 = t205 * t288;
t265 = t173 * t281;
t263 = qJD(1) * t292;
t195 = -t202 * t228 - t227 * t294;
t247 = t202 * t227 - t228 * t294;
t166 = 0.1e1 / t168;
t161 = t299 * t204;
t155 = t248 * t268 + t311;
t154 = 0.2e1 * t248 * t288 - t311;
t1 = [t266 * t295 + (t176 * t205 + t204 * t303) * t184, t154, 0, 0, t155, 0; -t200 * t169 * t271 + (-t177 * t169 + (-t157 * t200 - t161 * t176) * t170) * t166 - ((-t161 * t270 + t299 * t282) * t166 + (-t161 * t271 + ((t162 * t184 * t276 + t268) * t280 + (0.2e1 * t200 * t266 - t162 + (t162 - t250) * t184) * t279) * t166) * t170) * t204, t319 + (((t154 * t200 + t314) * t279 + (-t154 * t208 + t313) * t280 - t317) * t170 - t316) * t166, 0, 0, -t319 + (((t155 * t200 - t314) * t279 + (-t155 * t208 - t313) * t280 + t317) * t170 + t316) * t166, 0; (t267 * t296 - t265) * t195 - (-t174 * t281 + t191 * t269) * t247 + ((qJD(6) * t195 - t178 * t227 - t228 * t263) * t191 + (qJD(6) * t247 - t178 * t228 + t227 * t263) * t277 + t195 * t260) * t179, t249 * t289 * t295 + t318, 0, 0, -t249 * t204 * t269 - t318, t269 + (t265 - (-t179 * t284 - t267) * t246) * t296;];
JaD_rot  = t1;
