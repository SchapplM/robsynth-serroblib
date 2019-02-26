% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR9
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR9_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR9_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:12
% EndTime: 2019-02-26 22:08:12
% DurationCPUTime: 0.30s
% Computational Cost: add. (344->68), mult. (1036->114), div. (0->0), fcn. (1016->10), ass. (0->47)
t325 = cos(pkin(6));
t341 = qJD(2) * t325 + qJD(1);
t327 = sin(qJ(2));
t328 = sin(qJ(1));
t355 = t328 * t327;
t348 = t325 * t355;
t330 = cos(qJ(2));
t331 = cos(qJ(1));
t352 = t331 * t330;
t307 = -qJD(1) * t348 - qJD(2) * t355 + t341 * t352;
t323 = sin(pkin(6));
t356 = t323 * t331;
t363 = -qJD(3) * t356 + t307;
t353 = t331 * t327;
t354 = t328 * t330;
t311 = t325 * t353 + t354;
t351 = qJD(1) * t323;
t362 = -qJD(3) * t311 + t328 * t351;
t326 = sin(qJ(3));
t329 = cos(qJ(3));
t361 = t362 * t326 + t363 * t329;
t322 = sin(pkin(11));
t324 = cos(pkin(11));
t340 = t324 * r_i_i_C(1) - t322 * r_i_i_C(2) + pkin(3);
t359 = r_i_i_C(3) + qJ(4);
t360 = t359 * t326 + t340 * t329 + pkin(2);
t358 = t323 * t328;
t357 = t323 * t329;
t350 = qJD(2) * t330;
t347 = t325 * t352;
t345 = t331 * t351;
t344 = t323 * t350;
t339 = t322 * r_i_i_C(1) + t324 * r_i_i_C(2) + pkin(9);
t313 = -t348 + t352;
t338 = -t313 * t326 + t328 * t357;
t337 = t313 * t329 + t326 * t358;
t336 = t325 * t326 + t327 * t357;
t335 = t325 * t354 + t353;
t300 = t363 * t326 - t362 * t329;
t332 = t326 * qJD(4) + (-t340 * t326 + t359 * t329) * qJD(3);
t308 = t336 * qJD(3) + t326 * t344;
t306 = t335 * qJD(1) + t311 * qJD(2);
t305 = t311 * qJD(1) + t335 * qJD(2);
t304 = -qJD(1) * t347 - t331 * t350 + t341 * t355;
t299 = t338 * qJD(3) - t305 * t329 + t326 * t345;
t298 = t337 * qJD(3) - t305 * t326 - t329 * t345;
t1 = [(-t306 * t322 - t324 * t361) * r_i_i_C(1) + (-t306 * t324 + t322 * t361) * r_i_i_C(2) - t361 * pkin(3) - (t311 * t326 + t329 * t356) * qJD(4) - t307 * pkin(2) - t306 * pkin(9) - t359 * t300 + (-t331 * pkin(1) - pkin(8) * t358) * qJD(1), t360 * t304 - t339 * t305 - t332 * t335, t337 * qJD(4) - t340 * t298 + t359 * t299, t298, 0, 0; (t299 * t324 - t304 * t322) * r_i_i_C(1) + (-t299 * t322 - t304 * t324) * r_i_i_C(2) + t299 * pkin(3) - t338 * qJD(4) - t305 * pkin(2) - t304 * pkin(9) + t359 * t298 + (-t328 * pkin(1) + pkin(8) * t356) * qJD(1), t339 * t307 + t332 * (t347 - t355) - t360 * t306 -(-t311 * t329 + t326 * t356) * qJD(4) + t359 * t361 - t340 * t300, t300, 0, 0; 0 (t332 * t330 + (-t327 * t360 + t339 * t330) * qJD(2)) * t323, t336 * qJD(4) + t359 * (t329 * t344 + (-t323 * t326 * t327 + t325 * t329) * qJD(3)) - t340 * t308, t308, 0, 0;];
JaD_transl  = t1;
