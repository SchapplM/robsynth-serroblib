% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:25
% EndTime: 2019-02-26 22:28:25
% DurationCPUTime: 0.37s
% Computational Cost: add. (523->72), mult. (738->107), div. (0->0), fcn. (620->8), ass. (0->55)
t298 = qJ(3) + qJ(4);
t296 = cos(t298);
t345 = r_i_i_C(3) + qJ(5);
t357 = t345 * t296;
t295 = sin(t298);
t297 = qJD(3) + qJD(4);
t299 = sin(qJ(3));
t344 = pkin(3) * qJD(3);
t346 = r_i_i_C(1) + pkin(9) + pkin(8);
t308 = t346 * qJD(2) + qJD(5) * t295 - t299 * t344;
t348 = pkin(4) - r_i_i_C(2);
t356 = t308 - (t348 * t295 - t357) * t297;
t302 = cos(qJ(3));
t294 = t302 * pkin(3) + pkin(2);
t303 = cos(qJ(2));
t300 = sin(qJ(2));
t335 = qJD(2) * t300;
t355 = -t294 * t335 + t308 * t303;
t343 = t297 * t300;
t331 = t296 * t343;
t334 = qJD(2) * t303;
t353 = t295 * t334 + t331;
t321 = t345 * t295;
t352 = t348 * t296 + t294 + t321;
t301 = sin(qJ(1));
t304 = cos(qJ(1));
t339 = t304 * t296;
t315 = t301 * t295 + t303 * t339;
t341 = t301 * t303;
t351 = t295 * t341 + t339;
t347 = pkin(3) * t299;
t340 = t304 * t295;
t338 = qJD(1) * t301;
t337 = qJD(1) * t303;
t336 = qJD(1) * t304;
t333 = t296 * qJD(5);
t332 = t302 * t344;
t329 = t297 * t340;
t327 = pkin(7) + t347;
t323 = t301 * t335;
t322 = t304 * t335;
t320 = -qJD(3) + t337;
t319 = t353 * r_i_i_C(2) + t300 * t333 + t334 * t357;
t318 = (-qJD(3) * t303 + qJD(1)) * t302;
t316 = t332 - t333;
t314 = -t294 * t303 - t346 * t300 - pkin(1);
t313 = t301 * t297 * t296 + t295 * t336;
t273 = qJD(1) * t351 + t295 * t322 - t297 * t315;
t274 = t303 * t329 + (t301 * t337 + t322) * t296 - t313;
t312 = t315 * qJD(5) + t348 * t273 - t345 * t274;
t275 = -t295 * t323 - t296 * t338 + t313 * t303 - t329;
t276 = t315 * qJD(1) - t296 * t323 - t297 * t351;
t311 = -(-t296 * t341 + t340) * qJD(5) + t345 * t276 - t348 * t275;
t309 = qJD(2) * t352;
t1 = [t316 * t304 - t348 * t276 - t345 * t275 - t355 * t301 + (-t327 * t301 + t314 * t304) * qJD(1) (-t304 * t309 - t346 * t338) * t303 + (-t304 * t356 + t352 * t338) * t300 (t304 * t318 + (t320 * t301 + t322) * t299) * pkin(3) + t312, t312, -t273, 0; t316 * t301 - t348 * t274 - t345 * t273 + t355 * t304 + (t314 * t301 + t327 * t304) * qJD(1) (-t301 * t309 + t346 * t336) * t303 + (-t301 * t356 - t336 * t352) * t300 (t301 * t318 + (-t320 * t304 + t323) * t299) * pkin(3) + t311, t311, t275, 0; 0, -t300 * t309 + t356 * t303 (-pkin(4) * t295 - t347) * t334 + (-t332 + (-pkin(4) * t296 - t321) * t297) * t300 + t319, -pkin(4) * t331 + (-pkin(4) * t334 - t345 * t343) * t295 + t319, t353, 0;];
JaD_transl  = t1;
