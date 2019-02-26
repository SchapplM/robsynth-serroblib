% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:10
% EndTime: 2019-02-26 22:13:11
% DurationCPUTime: 0.66s
% Computational Cost: add. (495->104), mult. (1488->157), div. (0->0), fcn. (1380->8), ass. (0->66)
t299 = sin(qJ(2));
t298 = sin(qJ(3));
t302 = cos(qJ(3));
t301 = cos(qJ(5));
t297 = sin(qJ(5));
t350 = pkin(5) + r_i_i_C(1);
t323 = -t350 * t297 - qJ(4);
t312 = -t301 * r_i_i_C(2) + t323;
t347 = -t301 * pkin(5) - pkin(3) - pkin(4);
t324 = t297 * r_i_i_C(2) + t347;
t314 = t301 * r_i_i_C(1) - t324;
t307 = t312 * t298 - t314 * t302 - pkin(2);
t303 = cos(qJ(2));
t333 = -r_i_i_C(3) - qJ(6) - pkin(9) + pkin(8);
t354 = t333 * t303;
t366 = t307 * t299 + t354;
t315 = t297 * t298 + t301 * t302;
t316 = t297 * t302 - t298 * t301;
t365 = (t314 * t298 + t312 * t302) * qJD(3) - t333 * qJD(2) - t298 * qJD(4) + (t315 * r_i_i_C(2) + t350 * t316) * qJD(5);
t362 = (-t299 * pkin(2) + t354) * qJD(2) - t299 * qJD(6);
t304 = cos(qJ(1));
t343 = t304 * t302;
t300 = sin(qJ(1));
t345 = t300 * t303;
t279 = t298 * t345 + t343;
t344 = t304 * t298;
t280 = t302 * t345 - t344;
t320 = t279 * t301 - t280 * t297;
t357 = t320 * r_i_i_C(2);
t339 = qJD(2) * t303;
t355 = t316 * t339;
t351 = (qJD(3) - qJD(5)) * t299;
t337 = qJD(3) * t302;
t341 = qJD(1) * t304;
t310 = t298 * t341 + t300 * t337;
t336 = qJD(3) * t304;
t329 = t298 * t336;
t340 = qJD(2) * t299;
t332 = t300 * t340;
t342 = qJD(1) * t300;
t277 = -t298 * t332 - t302 * t342 + t310 * t303 - t329;
t274 = t277 * t301;
t346 = t300 * t298;
t338 = qJD(2) * t304;
t331 = qJD(3) * t346;
t330 = t299 * t338;
t328 = t302 * t336;
t327 = t297 * pkin(5) + qJ(4);
t313 = t315 * qJD(5);
t322 = -(t355 + (-t315 * qJD(3) + t313) * t299) * r_i_i_C(1) - (t315 * t339 + t316 * t351) * r_i_i_C(2);
t282 = t303 * t343 + t346;
t278 = t282 * qJD(1) - t302 * t332 - t303 * t331 - t328;
t321 = -t278 * t297 + t274;
t319 = t279 * t297 + t280 * t301;
t281 = -t300 * t302 + t303 * t344;
t318 = t281 * t301 - t282 * t297;
t317 = t281 * t297 + t282 * t301;
t309 = -pkin(2) * t303 - t333 * t299 - pkin(1);
t275 = t279 * qJD(1) + t298 * t330 - t303 * t328 - t331;
t308 = t318 * qJD(5) - t275 * t297;
t276 = t303 * t329 + (t303 * t342 + t330) * t302 - t310;
t270 = -t317 * qJD(5) - t275 * t301 + t276 * t297;
t306 = t307 * qJD(2) - qJD(6);
t305 = t365 * t299 + t306 * t303;
t271 = -t276 * t301 + t308;
t1 = [-t274 * r_i_i_C(2) - t279 * qJD(4) - t314 * t278 + t323 * t277 + (t319 * r_i_i_C(2) - t350 * t320) * qJD(5) - t362 * t300 + (-t300 * pkin(7) + t309 * t304) * qJD(1), t305 * t304 - t366 * t342, t282 * qJD(4) + t312 * t276 + t314 * t275 + (t318 * r_i_i_C(2) + t350 * t317) * qJD(5), -t275, -t271 * r_i_i_C(2) + t350 * t270, t299 * t342 - t303 * t338; t271 * r_i_i_C(1) + t270 * r_i_i_C(2) - t275 * qJ(4) + t281 * qJD(4) + t347 * t276 + t308 * pkin(5) + t362 * t304 + (pkin(7) * t304 + t309 * t300) * qJD(1), t305 * t300 + t366 * t341, -t274 * r_i_i_C(1) + t280 * qJD(4) - t312 * t278 + t324 * t277 + (t350 * t319 + t357) * qJD(5), t277, t321 * r_i_i_C(1) + (-t277 * t297 - t278 * t301) * r_i_i_C(2) + (-t319 * r_i_i_C(1) - t357) * qJD(5) + (-t319 * qJD(5) + t321) * pkin(5), -t299 * t341 - t300 * t339; 0, t306 * t299 - t365 * t303 (t347 * t298 + t327 * t302) * t339 + (qJD(4) * t302 + pkin(5) * t313 + (-t327 * t298 + t347 * t302) * qJD(3)) * t299 - t322, t298 * t339 + t299 * t337 (t315 * t351 - t355) * pkin(5) + t322, -t340;];
JaD_transl  = t1;
