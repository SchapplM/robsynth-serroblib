% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:37
% EndTime: 2019-02-26 22:11:38
% DurationCPUTime: 0.36s
% Computational Cost: add. (796->72), mult. (798->101), div. (0->0), fcn. (668->10), ass. (0->56)
t310 = qJ(3) + pkin(10);
t306 = qJ(5) + t310;
t303 = cos(t306);
t353 = r_i_i_C(3) + qJ(6);
t366 = t353 * t303;
t309 = qJD(3) + qJD(5);
t296 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t310);
t285 = t296 * qJD(3);
t302 = sin(t306);
t354 = r_i_i_C(2) + pkin(9) + qJ(4) + pkin(8);
t324 = qJD(2) * t354 + qJD(6) * t302 - t285;
t356 = pkin(5) + r_i_i_C(1);
t340 = t356 * t302;
t365 = (t340 - t366) * t309 - t324;
t297 = pkin(4) * cos(t310) + cos(qJ(3)) * pkin(3);
t293 = pkin(2) + t297;
t312 = sin(qJ(2));
t315 = cos(qJ(2));
t364 = t324 * t315 + (pkin(7) + t296) * qJD(1) - (qJD(2) * t293 - qJD(4)) * t312;
t326 = -t302 * t353 - t303 * t356;
t321 = -t293 + t326;
t362 = t312 * t321 + t315 * t354;
t313 = sin(qJ(1));
t316 = cos(qJ(1));
t348 = t316 * t303;
t328 = t302 * t313 + t315 * t348;
t350 = t313 * t315;
t358 = t302 * t350 + t348;
t352 = t309 * t312;
t349 = t316 * t302;
t347 = qJD(1) * t313;
t346 = qJD(1) * t315;
t345 = qJD(1) * t316;
t344 = qJD(2) * t312;
t343 = qJD(2) * t315;
t342 = qJD(2) * t316;
t341 = t303 * qJD(6);
t339 = t303 * t352;
t337 = t309 * t349;
t335 = t312 * t341 + t343 * t366;
t333 = t313 * t344;
t332 = t312 * t342;
t329 = t296 * t346 - t285;
t327 = t303 * t309 * t313 + t302 * t345;
t277 = qJD(1) * t358 + t302 * t332 - t309 * t328;
t278 = t315 * t337 + (t313 * t346 + t332) * t303 - t327;
t323 = t328 * qJD(6) + t277 * t356 - t278 * t353;
t279 = -t302 * t333 - t303 * t347 + t315 * t327 - t337;
t280 = t328 * qJD(1) - t303 * t333 - t309 * t358;
t322 = -(-t303 * t350 + t349) * qJD(6) + t353 * t280 - t356 * t279;
t286 = t297 * qJD(3);
t320 = qJD(1) * t297 - t286 * t315 + t296 * t344;
t319 = qJD(2) * t321 + qJD(4);
t318 = -t341 + t286 + (-t293 * t315 - t312 * t354 - pkin(1)) * qJD(1);
t317 = t312 * t365 + t315 * t319;
t1 = [-t279 * t353 - t280 * t356 - t313 * t364 + t316 * t318, t316 * t317 - t347 * t362, t313 * t329 + t316 * t320 + t323, -t312 * t347 + t315 * t342, t323, -t277; -t277 * t353 - t278 * t356 + t313 * t318 + t316 * t364, t313 * t317 + t345 * t362, t313 * t320 - t316 * t329 + t322, t312 * t345 + t313 * t343, t322, t279; 0, t312 * t319 - t315 * t365 (-t296 - t340) * t343 + (t309 * t326 - t286) * t312 + t335, t344, -t356 * t339 + (-t343 * t356 - t352 * t353) * t302 + t335, t302 * t343 + t339;];
JaD_transl  = t1;
