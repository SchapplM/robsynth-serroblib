% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:08
% EndTime: 2019-02-26 20:11:08
% DurationCPUTime: 0.27s
% Computational Cost: add. (471->52), mult. (787->91), div. (0->0), fcn. (763->12), ass. (0->47)
t362 = r_i_i_C(3) + qJ(5);
t327 = sin(pkin(12));
t330 = cos(pkin(12));
t347 = r_i_i_C(1) * t330 - r_i_i_C(2) * t327 + pkin(4);
t328 = sin(pkin(11));
t331 = cos(pkin(11));
t332 = cos(pkin(6));
t336 = cos(qJ(2));
t354 = t332 * t336;
t350 = t331 * t354;
t334 = sin(qJ(2));
t353 = qJD(2) * t334;
t310 = -qJD(2) * t350 + t328 * t353;
t325 = qJD(3) + qJD(4);
t329 = sin(pkin(6));
t358 = t329 * t331;
t366 = -t325 * t358 - t310;
t326 = qJ(3) + qJ(4);
t323 = sin(t326);
t324 = cos(t326);
t335 = cos(qJ(3));
t365 = pkin(3) * t335 + t362 * t323 + t347 * t324 + pkin(2);
t361 = t323 * t325;
t360 = t324 * t325;
t359 = t328 * t329;
t333 = sin(qJ(3));
t357 = t329 * t333;
t356 = t329 * t334;
t355 = t332 * t334;
t351 = t325 * t356;
t349 = qJD(2) * t329 * t336;
t345 = t328 * t354 + t331 * t334;
t312 = t345 * qJD(2);
t348 = t325 * t359 - t312;
t346 = t327 * r_i_i_C(1) + t330 * r_i_i_C(2) + pkin(8) + pkin(9);
t315 = t328 * t336 + t331 * t355;
t344 = t328 * t355 - t331 * t336;
t343 = t325 * t332 + t349;
t296 = t315 * t360 + t366 * t323;
t342 = -qJD(5) * (-t315 * t324 + t323 * t358) + t362 * (-t315 * t361 + t366 * t324) - t347 * t296;
t298 = t348 * t323 - t344 * t360;
t341 = -qJD(5) * (-t323 * t359 + t324 * t344) + t362 * (t348 * t324 + t344 * t361) - t347 * t298;
t303 = t343 * t323 + t324 * t351;
t340 = -(-t323 * t332 - t324 * t356) * qJD(5) + t362 * (-t323 * t351 + t343 * t324) - t347 * t303;
t339 = qJD(2) * t365;
t338 = -pkin(3) * qJD(3) * t333 + qJD(5) * t323 + (-t347 * t323 + t362 * t324) * t325;
t1 = [0, -t346 * t312 - t338 * t345 + t344 * t339 (t312 * t333 + (-t328 * t357 + t335 * t344) * qJD(3)) * pkin(3) + t341, t341, t298, 0; 0, -t346 * t310 - t315 * t339 + t338 * (-t328 * t334 + t350) (t310 * t333 + (-t315 * t335 + t331 * t357) * qJD(3)) * pkin(3) + t342, t342, t296, 0; 0 (-t365 * t353 + (t346 * qJD(2) + t338) * t336) * t329 (-t333 * t349 + (-t332 * t333 - t335 * t356) * qJD(3)) * pkin(3) + t340, t340, t303, 0;];
JaD_transl  = t1;
