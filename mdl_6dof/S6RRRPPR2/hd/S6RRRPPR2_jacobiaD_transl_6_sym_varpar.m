% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:57
% EndTime: 2019-02-26 22:03:57
% DurationCPUTime: 0.30s
% Computational Cost: add. (503->69), mult. (493->90), div. (0->0), fcn. (370->10), ass. (0->60)
t271 = qJ(2) + qJ(3);
t266 = pkin(10) + t271;
t265 = cos(t266);
t272 = sin(qJ(6));
t275 = cos(qJ(6));
t326 = r_i_i_C(1) * t272 + r_i_i_C(2) * t275 + qJ(5);
t284 = t326 * t265;
t264 = sin(t266);
t273 = sin(qJ(2));
t315 = pkin(2) * qJD(2);
t304 = t273 * t315;
t267 = sin(t271);
t270 = qJD(2) + qJD(3);
t314 = t265 * t270;
t319 = pkin(3) * t270;
t323 = -qJ(5) * t314 + t267 * t319;
t306 = pkin(4) + pkin(9) + r_i_i_C(3);
t327 = -t306 * t270 + qJD(5);
t329 = t327 * t264 + (pkin(5) + qJ(4) + pkin(8) + pkin(7)) * qJD(1) - t304 - t323;
t307 = qJD(6) * t275;
t297 = t265 * t307;
t328 = r_i_i_C(1) * t297 + qJD(5) * t265;
t291 = qJD(6) * t264 + qJD(1);
t325 = t275 * t291;
t324 = t291 * t272;
t321 = pkin(3) * t267;
t268 = cos(t271);
t320 = pkin(3) * t268;
t313 = t270 * t272;
t312 = t270 * t275;
t274 = sin(qJ(1));
t311 = qJD(1) * t274;
t277 = cos(qJ(1));
t310 = qJD(1) * t277;
t308 = qJD(6) * t272;
t303 = t265 * t312;
t302 = t274 * t314;
t301 = t277 * t314;
t299 = t264 * t311;
t298 = t265 * t308;
t296 = t306 * t264;
t295 = t306 * t265;
t292 = r_i_i_C(2) * t298;
t290 = -qJD(1) * t264 - qJD(6);
t288 = t328 * t274 + t310 * t284;
t287 = t328 * t277 + t306 * t299;
t286 = t290 * t277;
t283 = t290 * t274 + t301;
t276 = cos(qJ(2));
t282 = qJD(4) + (-t276 * pkin(2) - qJ(5) * t264 - pkin(1) - t295 - t320) * qJD(1);
t281 = -t264 * t326 - t295;
t280 = t265 * r_i_i_C(1) * t313 + r_i_i_C(2) * t303 - t323 + (r_i_i_C(1) * t307 - r_i_i_C(2) * t308 + t327) * t264;
t279 = -t268 * t319 + t281 * t270 - t276 * t315 - t292;
t278 = -t292 + (t281 - t320) * t270;
t260 = -t273 * pkin(2) - t321;
t244 = t283 * t272 + t277 * t325;
t243 = t283 * t275 - t277 * t324;
t242 = -t274 * t325 + (t286 - t302) * t272;
t241 = t275 * t286 + (-t303 + t324) * t274;
t1 = [t242 * r_i_i_C(1) + t241 * r_i_i_C(2) - t329 * t274 + t282 * t277 (-t260 - t284) * t311 + t279 * t277 + t287 (-t284 + t321) * t311 + t278 * t277 + t287, t310, -t299 + t301, t243 * r_i_i_C(1) - t244 * r_i_i_C(2); t244 * r_i_i_C(1) + t243 * r_i_i_C(2) + t282 * t274 + t329 * t277 (t260 - t296) * t310 + t279 * t274 + t288 (-t296 - t321) * t310 + t278 * t274 + t288, t311, t264 * t310 + t302, -t241 * r_i_i_C(1) + t242 * r_i_i_C(2); 0, t280 - t304, t280, 0, t270 * t264 (-t264 * t313 + t297) * r_i_i_C(2) + (t264 * t312 + t298) * r_i_i_C(1);];
JaD_transl  = t1;
