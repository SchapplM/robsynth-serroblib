% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiaD_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:06
% EndTime: 2019-02-26 19:44:07
% DurationCPUTime: 0.29s
% Computational Cost: add. (257->72), mult. (859->144), div. (0->0), fcn. (972->14), ass. (0->56)
t318 = sin(pkin(14));
t322 = sin(pkin(6));
t323 = cos(pkin(14));
t329 = sin(qJ(3));
t331 = cos(qJ(3));
t326 = cos(pkin(7));
t343 = t326 * t329;
t321 = sin(pkin(7));
t327 = cos(pkin(6));
t349 = t321 * t327;
t304 = (t318 * t331 + t323 * t343) * t322 + t329 * t349;
t324 = cos(pkin(13));
t319 = sin(pkin(13));
t353 = t319 * t327;
t313 = -t318 * t353 + t323 * t324;
t312 = -t318 * t324 - t323 * t353;
t350 = t321 * t322;
t334 = t312 * t326 + t319 * t350;
t300 = t313 * t331 + t334 * t329;
t346 = t324 * t327;
t311 = t318 * t346 + t319 * t323;
t356 = t311 * t329;
t355 = t311 * t331;
t320 = sin(pkin(8));
t328 = sin(qJ(4));
t352 = t320 * t328;
t330 = cos(qJ(4));
t351 = t320 * t330;
t348 = t322 * t323;
t347 = t322 * t326;
t325 = cos(pkin(8));
t345 = t325 * t328;
t344 = t325 * t330;
t342 = qJD(3) * t329;
t341 = qJD(3) * t331;
t340 = t324 * t350;
t338 = t321 * t341;
t337 = t326 * t341;
t336 = r_i_i_C(1) * t330 - r_i_i_C(2) * t328 + pkin(3);
t310 = -t318 * t319 + t323 * t346;
t335 = -t310 * t326 + t340;
t332 = (r_i_i_C(1) * t328 + r_i_i_C(2) * t330) * t325 + (-pkin(10) - r_i_i_C(3)) * t320;
t309 = -t321 * t348 + t326 * t327;
t306 = -t312 * t321 + t319 * t347;
t305 = -t310 * t321 - t324 * t347;
t303 = t331 * t349 + (t323 * t326 * t331 - t318 * t329) * t322;
t302 = t304 * qJD(3);
t301 = t318 * t322 * t342 - t327 * t338 - t337 * t348;
t299 = -t313 * t329 + t334 * t331;
t298 = t310 * t343 - t329 * t340 + t355;
t297 = -t335 * t331 - t356;
t296 = t300 * qJD(3);
t295 = -t319 * t322 * t338 - t312 * t337 + t313 * t342;
t294 = (t335 * t329 - t355) * qJD(3);
t293 = -t310 * t337 + (t331 * t340 + t356) * qJD(3);
t1 = [0, 0, -t336 * t296 + t332 * t295 + ((-t299 * t328 - t300 * t344) * r_i_i_C(1) + (-t299 * t330 + t300 * t345) * r_i_i_C(2)) * qJD(4) (t295 * t328 - t296 * t344) * r_i_i_C(1) + (t295 * t330 + t296 * t345) * r_i_i_C(2) + ((-t299 * t345 - t300 * t330 - t306 * t352) * r_i_i_C(1) + (-t299 * t344 + t300 * t328 - t306 * t351) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, t336 * t294 + t332 * t293 + ((-t297 * t328 - t298 * t344) * r_i_i_C(1) + (-t297 * t330 + t298 * t345) * r_i_i_C(2)) * qJD(4) (t293 * t328 + t294 * t344) * r_i_i_C(1) + (t293 * t330 - t294 * t345) * r_i_i_C(2) + ((-t297 * t345 - t298 * t330 - t305 * t352) * r_i_i_C(1) + (-t297 * t344 + t298 * t328 - t305 * t351) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, -t336 * t302 + t332 * t301 + ((-t303 * t328 - t304 * t344) * r_i_i_C(1) + (-t303 * t330 + t304 * t345) * r_i_i_C(2)) * qJD(4) (t301 * t328 - t302 * t344) * r_i_i_C(1) + (t301 * t330 + t302 * t345) * r_i_i_C(2) + ((-t303 * t345 - t304 * t330 - t309 * t352) * r_i_i_C(1) + (-t303 * t344 + t304 * t328 - t309 * t351) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
