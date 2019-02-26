% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP3
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
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:37
% EndTime: 2019-02-26 22:10:37
% DurationCPUTime: 0.44s
% Computational Cost: add. (709->85), mult. (745->119), div. (0->0), fcn. (606->10), ass. (0->62)
t319 = qJ(2) + qJ(3);
t315 = sin(t319);
t373 = pkin(5) + r_i_i_C(1);
t381 = t315 * t373;
t317 = pkin(10) + qJ(5);
t314 = cos(t317);
t313 = sin(t317);
t325 = cos(qJ(1));
t357 = qJD(5) * t325;
t346 = t313 * t357;
t323 = sin(qJ(1));
t362 = qJD(1) * t323;
t380 = t314 * t362 + t346;
t311 = cos(pkin(10)) * pkin(4) + pkin(3);
t369 = r_i_i_C(3) + qJ(6);
t379 = t369 * t313 + t311;
t316 = cos(t319);
t350 = t316 * t362;
t318 = qJD(2) + qJD(3);
t363 = t318 * t325;
t353 = t315 * t363;
t377 = t350 + t353;
t322 = sin(qJ(2));
t368 = pkin(2) * qJD(2);
t355 = t322 * t368;
t356 = qJD(6) * t313;
t321 = -pkin(9) - qJ(4);
t370 = r_i_i_C(2) - t321;
t376 = (t370 * t318 + t356) * t316 + (pkin(4) * sin(pkin(10)) + pkin(8) + pkin(7)) * qJD(1) - (t311 * t318 - qJD(4)) * t315 - t355;
t372 = pkin(2) * t322;
t371 = r_i_i_C(2) * t316;
t367 = t314 * t325;
t366 = t316 * t318;
t365 = t316 * t323;
t364 = t318 * t323;
t361 = qJD(1) * t325;
t360 = qJD(4) * t316;
t359 = qJD(5) * t314;
t358 = qJD(5) * t323;
t354 = t315 * t364;
t352 = t373 * qJD(5);
t351 = t315 * t362;
t347 = t313 * t358;
t345 = t314 * t357;
t343 = qJD(5) * t369;
t337 = t321 * t354 + t323 * t360 + t347 * t381 + t361 * t371;
t336 = t313 * t323 + t316 * t367;
t335 = t313 * t361 + t314 * t358;
t334 = -t314 * t343 - t356;
t333 = -t373 * t314 - t379;
t332 = t377 * t321 + t325 * t360 + t379 * t351 + t380 * t381;
t324 = cos(qJ(2));
t331 = -t314 * qJD(6) + (-pkin(2) * t324 - t311 * t316 - t370 * t315 - pkin(1)) * qJD(1);
t330 = t333 * t315 - t316 * t321;
t329 = t334 * t315 + (-r_i_i_C(2) * t315 + t333 * t316) * t318;
t328 = r_i_i_C(2) * t366 + t315 * qJD(4) + t330 * t318 + (-t313 * t352 + t369 * t359 + t356) * t316;
t327 = -t324 * t368 + t329;
t280 = t336 * qJD(1) - t314 * t354 - t316 * t347 - t345;
t279 = -t313 * t354 + t335 * t316 - t380;
t278 = t377 * t314 + t316 * t346 - t335;
t277 = t313 * t353 - t316 * t345 - t347 + (t313 * t365 + t367) * qJD(1);
t1 = [-t369 * t279 - t373 * t280 - t376 * t323 + t331 * t325 (-t371 + t372) * t362 + t327 * t325 + t332, -r_i_i_C(2) * t350 + t329 * t325 + t332, t316 * t363 - t351, t336 * qJD(6) + t373 * t277 - t369 * t278, -t277; -t369 * t277 - t373 * t278 + t331 * t323 + t376 * t325 (t330 - t372) * t361 + t327 * t323 + t337 (-t321 * t361 + t333 * t364) * t316 + ((-r_i_i_C(2) * t318 + t334) * t323 + t333 * t361) * t315 + t337, t315 * t361 + t316 * t364 -(t313 * t325 - t314 * t365) * qJD(6) + t369 * t280 - t373 * t279, t279; 0, t328 - t355, t328, t318 * t315 (-t315 * t343 - t373 * t366) * t313 + (t369 * t366 + (qJD(6) - t352) * t315) * t314, t313 * t366 + t315 * t359;];
JaD_transl  = t1;
