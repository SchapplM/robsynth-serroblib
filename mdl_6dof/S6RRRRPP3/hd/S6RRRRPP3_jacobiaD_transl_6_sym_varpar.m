% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP3
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
% Datum: 2019-02-26 22:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:23
% EndTime: 2019-02-26 22:26:23
% DurationCPUTime: 0.46s
% Computational Cost: add. (635->81), mult. (928->113), div. (0->0), fcn. (759->8), ass. (0->61)
t322 = sin(qJ(4));
t375 = r_i_i_C(2) + qJ(5);
t381 = t375 * t322;
t325 = cos(qJ(4));
t327 = cos(qJ(1));
t362 = qJD(4) * t327;
t324 = sin(qJ(1));
t366 = qJD(1) * t324;
t339 = t322 * t362 + t325 * t366;
t321 = qJ(2) + qJ(3);
t318 = sin(t321);
t320 = qJD(2) + qJD(3);
t319 = cos(t321);
t360 = pkin(5) + pkin(9) + r_i_i_C(1);
t347 = t360 * t319;
t323 = sin(qJ(2));
t374 = pkin(2) * qJD(2);
t358 = t323 * t374;
t380 = (-pkin(3) * t318 + t347) * t320 - t358;
t359 = pkin(4) + r_i_i_C(3) + qJ(6);
t378 = t318 * t359;
t377 = pkin(2) * t323;
t373 = t318 * t320;
t372 = t319 * t320;
t371 = t319 * t325;
t370 = t324 * t322;
t369 = t324 * t325;
t368 = t327 * t322;
t367 = t327 * t325;
t365 = qJD(1) * t327;
t364 = qJD(4) * t322;
t363 = qJD(4) * t325;
t361 = qJD(5) * t322;
t357 = t324 * t373;
t356 = t327 * t373;
t351 = t324 * t364;
t349 = t325 * t362;
t348 = t360 * t318;
t346 = t359 * t322;
t341 = t365 * t347 + t351 * t378;
t340 = t319 * t367 + t370;
t281 = t319 * t370 + t367;
t338 = t322 * t365 + t324 * t363;
t326 = cos(qJ(2));
t337 = -t326 * pkin(2) - pkin(3) * t319 - pkin(1) - t348;
t336 = -t359 * t325 - t381;
t335 = t339 * t378 + (pkin(3) + t381) * t318 * t366;
t334 = -pkin(3) + t336;
t333 = -t361 + (-t375 * qJD(4) - qJD(6)) * t325;
t332 = t334 * t320;
t331 = qJD(6) * t371 + t318 * t332 + t360 * t372 + (-qJD(4) * t346 + t375 * t363 + t361) * t319;
t330 = t333 * t318 + (t334 * t319 - t348) * t320;
t329 = -t326 * t374 + t330;
t328 = -pkin(8) - pkin(7);
t283 = t319 * t368 - t369;
t282 = t319 * t369 - t368;
t280 = t340 * qJD(1) - t319 * t351 - t325 * t357 - t349;
t279 = t338 * t319 - t322 * t357 - t339;
t278 = t339 * t319 + t325 * t356 - t338;
t277 = t281 * qJD(1) - t319 * t349 + t322 * t356 - t351;
t1 = [-t281 * qJD(5) - t282 * qJD(6) - t375 * t279 - t359 * t280 - t380 * t324 + (t324 * t328 + t337 * t327) * qJD(1) (-t347 + t377) * t366 + t329 * t327 + t335, t330 * t327 - t347 * t366 + t335, qJD(5) * t340 - t283 * qJD(6) + t359 * t277 - t375 * t278, -t277, -t278; t283 * qJD(5) + t340 * qJD(6) - t375 * t277 - t359 * t278 + t380 * t327 + (t337 * t324 - t327 * t328) * qJD(1) (t334 * t318 - t377) * t365 + t329 * t324 + t341, t324 * t319 * t332 + (t334 * t365 + (-t360 * t320 + t333) * t324) * t318 + t341, t282 * qJD(5) - t281 * qJD(6) - t359 * t279 + t375 * t280, t279, t280; 0, t331 - t358, t331 (t375 * t325 - t346) * t372 + (t336 * qJD(4) + qJD(5) * t325 - qJD(6) * t322) * t318, t318 * t363 + t322 * t372, -t318 * t364 + t320 * t371;];
JaD_transl  = t1;
