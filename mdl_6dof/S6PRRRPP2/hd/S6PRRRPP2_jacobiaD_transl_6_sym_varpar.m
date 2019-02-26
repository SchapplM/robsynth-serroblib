% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:09:34
% EndTime: 2019-02-26 20:09:34
% DurationCPUTime: 0.49s
% Computational Cost: add. (560->95), mult. (1745->161), div. (0->0), fcn. (1809->10), ass. (0->63)
t328 = sin(qJ(3));
t331 = cos(qJ(3));
t359 = pkin(9) - r_i_i_C(3) - qJ(6);
t327 = sin(qJ(4));
t362 = qJD(5) * t327;
t377 = (-pkin(3) * t328 + t359 * t331) * qJD(3) - t328 * qJD(6) + t331 * t362;
t329 = sin(qJ(2));
t332 = cos(qJ(2));
t325 = sin(pkin(10));
t371 = cos(pkin(6));
t353 = t325 * t371;
t370 = cos(pkin(10));
t315 = -t329 * t353 + t370 * t332;
t326 = sin(pkin(6));
t367 = t326 * t331;
t317 = t371 * t328 + t329 * t367;
t330 = cos(qJ(4));
t366 = t326 * t332;
t376 = -t317 * t330 + t327 * t366;
t364 = qJD(3) * t328;
t375 = (qJD(2) * t331 - qJD(4)) * t329 + t332 * t364;
t372 = r_i_i_C(2) + qJ(5);
t368 = t326 * t328;
t365 = qJD(2) * t329;
t363 = qJD(4) * t331;
t361 = qJD(5) * t330;
t360 = pkin(4) + pkin(5) + r_i_i_C(1);
t357 = t329 * t368;
t356 = t326 * t365;
t355 = qJD(2) * t366;
t352 = t326 * t370;
t348 = t331 * t352;
t345 = t371 * t370;
t341 = t332 * t345;
t308 = -qJD(2) * t341 + t325 * t365;
t312 = t325 * t329 - t341;
t347 = t312 * t363 - t308;
t314 = t370 * t329 + t332 * t353;
t310 = t314 * qJD(2);
t346 = t314 * t363 - t310;
t313 = t325 * t332 + t329 * t345;
t339 = -t313 * t331 + t328 * t352;
t344 = t312 * t327 - t330 * t339;
t301 = t315 * t331 + t325 * t368;
t343 = t301 * t330 + t314 * t327;
t342 = (qJD(2) - t363) * t332;
t338 = -t331 * pkin(3) - t359 * t328 - pkin(2);
t337 = t372 * t327 + t360 * t330 + pkin(3);
t309 = t313 * qJD(2);
t336 = qJD(4) * t313 - t309 * t331 + t312 * t364;
t311 = t315 * qJD(2);
t335 = qJD(4) * t315 - t311 * t331 + t314 * t364;
t334 = t362 + (-t360 * t327 + t372 * t330) * qJD(4);
t303 = -qJD(3) * t357 + (t371 * qJD(3) + t355) * t331;
t302 = -t317 * qJD(3) - t328 * t355;
t297 = -t315 * t364 + (qJD(3) * t325 * t326 - t310) * t331;
t296 = -t301 * qJD(3) + t310 * t328;
t295 = -qJD(3) * t348 - t308 * t331 - t313 * t364;
t294 = t339 * qJD(3) + t308 * t328;
t290 = -t376 * qJD(4) + t303 * t327 - t330 * t356;
t284 = t343 * qJD(4) + t297 * t327 - t311 * t330;
t282 = t344 * qJD(4) + t295 * t327 - t309 * t330;
t1 = [0, -t315 * t361 - t310 * pkin(8) + t372 * (t335 * t327 - t346 * t330) + t360 * (t346 * t327 + t335 * t330) - t377 * t314 + t338 * t311, -t301 * qJD(6) + t359 * t297 + t334 * (-t315 * t328 + t325 * t367) + t337 * t296, t343 * qJD(5) + t372 * (t297 * t330 + t311 * t327 + (-t301 * t327 + t314 * t330) * qJD(4)) - t360 * t284, t284, t296; 0, -t313 * t361 - t308 * pkin(8) + t372 * (t336 * t327 - t347 * t330) + t360 * (t347 * t327 + t336 * t330) - t377 * t312 + t338 * t309, t339 * qJD(6) + t359 * t295 + t334 * (-t313 * t328 - t348) + t337 * t294, t344 * qJD(5) + t372 * (t295 * t330 + t309 * t327 + (t312 * t330 + t327 * t339) * qJD(4)) - t360 * t282, t282, t294; 0 (-t372 * (t375 * t327 + t330 * t342) + t360 * (t327 * t342 - t375 * t330) - t329 * t361 + t377 * t332 + (t332 * pkin(8) + t338 * t329) * qJD(2)) * t326, -t317 * qJD(6) + t359 * t303 + t334 * (t371 * t331 - t357) + t337 * t302, -t376 * qJD(5) + t372 * (t327 * t356 + t303 * t330 + (-t317 * t327 - t330 * t366) * qJD(4)) - t360 * t290, t290, t302;];
JaD_transl  = t1;
