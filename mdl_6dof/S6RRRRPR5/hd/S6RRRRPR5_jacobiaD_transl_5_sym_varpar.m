% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:54
% EndTime: 2019-02-26 22:32:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (474->71), mult. (702->99), div. (0->0), fcn. (565->8), ass. (0->56)
t316 = sin(qJ(4));
t319 = cos(qJ(4));
t321 = cos(qJ(1));
t357 = qJD(4) * t321;
t318 = sin(qJ(1));
t360 = qJD(1) * t318;
t330 = t316 * t357 + t319 * t360;
t369 = -r_i_i_C(1) - pkin(4);
t367 = r_i_i_C(3) + qJ(5);
t374 = t367 * t316;
t315 = qJ(2) + qJ(3);
t313 = cos(t315);
t370 = pkin(9) + r_i_i_C(2);
t352 = t370 * t313;
t314 = qJD(2) + qJD(3);
t351 = t370 * t314;
t317 = sin(qJ(2));
t366 = pkin(2) * qJD(2);
t354 = t317 * t366;
t356 = qJD(5) * t316;
t312 = sin(t315);
t365 = t312 * t314;
t373 = -pkin(3) * t365 + (t351 + t356) * t313 - t354;
t358 = qJD(4) * t319;
t328 = -t367 * t358 - t356;
t368 = pkin(2) * t317;
t364 = t313 * t314;
t363 = t313 * t318;
t362 = t318 * t316;
t361 = t321 * t319;
t359 = qJD(1) * t321;
t355 = t319 * qJD(5);
t353 = t370 * t312;
t350 = t369 * t316;
t349 = t318 * t365;
t348 = t321 * t365;
t343 = qJD(4) * t362;
t341 = t319 * t357;
t340 = -t369 * t312 * t343 + t359 * t352;
t335 = t313 * t361 + t362;
t320 = cos(qJ(2));
t333 = -t320 * pkin(2) - pkin(3) * t313 - pkin(1) - t353;
t332 = (-t369 * t330 + (pkin(3) + t374) * t360) * t312;
t331 = t369 * t319 - t374;
t329 = t316 * t359 + t318 * t358;
t327 = -pkin(3) + t331;
t326 = t327 * t314;
t325 = t312 * t326 + t370 * t364 + (qJD(4) * t350 - t328) * t313;
t324 = t328 * t312 + (t327 * t313 - t353) * t314;
t323 = -t320 * t366 + t324;
t322 = -pkin(8) - pkin(7);
t284 = t335 * qJD(1) - t313 * t343 - t319 * t349 - t341;
t283 = t329 * t313 - t316 * t349 - t330;
t282 = t330 * t313 + t319 * t348 - t329;
t281 = t316 * t348 - t313 * t341 - t343 + (t313 * t362 + t361) * qJD(1);
t1 = [-t321 * t355 + t369 * t284 - t367 * t283 - t373 * t318 + (t318 * t322 + t333 * t321) * qJD(1) (-t352 + t368) * t360 + t323 * t321 + t332, t324 * t321 - t352 * t360 + t332, t335 * qJD(5) - t369 * t281 - t367 * t282, -t281, 0; -t318 * t355 + t369 * t282 - t367 * t281 + t373 * t321 + (t333 * t318 - t321 * t322) * qJD(1) (t327 * t312 - t368) * t359 + t323 * t318 + t340, t326 * t363 + ((-t351 + t328) * t318 + t327 * t359) * t312 + t340 -(t321 * t316 - t319 * t363) * qJD(5) + t367 * t284 + t369 * t283, t283, 0; 0, t325 - t354, t325 (t367 * t319 + t350) * t364 + (t331 * qJD(4) + t355) * t312, t312 * t358 + t316 * t364, 0;];
JaD_transl  = t1;
