% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:02
% EndTime: 2019-02-26 21:43:02
% DurationCPUTime: 0.30s
% Computational Cost: add. (429->69), mult. (881->107), div. (0->0), fcn. (857->10), ass. (0->51)
t298 = pkin(11) + qJ(4);
t296 = sin(t298);
t297 = cos(t298);
t335 = r_i_i_C(3) + qJ(5);
t337 = -r_i_i_C(2) + pkin(4);
t342 = (t337 * t296 - t335 * t297) * qJD(4) - qJD(5) * t296;
t304 = sin(qJ(1));
t301 = cos(pkin(6));
t316 = qJD(2) * t301 + qJD(1);
t303 = sin(qJ(2));
t330 = t304 * t303;
t322 = t301 * t330;
t325 = qJD(2) * t303;
t305 = cos(qJ(2));
t306 = cos(qJ(1));
t327 = t306 * t305;
t281 = -qJD(1) * t322 - t304 * t325 + t316 * t327;
t328 = t306 * t303;
t329 = t304 * t305;
t286 = t301 * t328 + t329;
t300 = sin(pkin(6));
t331 = t300 * t306;
t315 = t286 * t296 + t297 * t331;
t326 = qJD(1) * t300;
t319 = t304 * t326;
t341 = t315 * qJD(4) - t281 * t297 - t296 * t319;
t314 = -t286 * t297 + t296 * t331;
t340 = t314 * qJD(4) - t281 * t296 + t297 * t319;
t295 = cos(pkin(11)) * pkin(3) + pkin(2);
t338 = t335 * t296 + t337 * t297 + t295;
t336 = r_i_i_C(1) + pkin(9) + qJ(3);
t288 = -t322 + t327;
t334 = t288 * t296;
t333 = t300 * t303;
t332 = t300 * t304;
t324 = qJD(2) * t305;
t321 = t301 * t327;
t320 = pkin(3) * sin(pkin(11)) + pkin(8);
t318 = t306 * t326;
t317 = t300 * t324;
t313 = t288 * t297 + t296 * t332;
t312 = t301 * t296 + t297 * t333;
t311 = t301 * t329 + t328;
t285 = -t321 + t330;
t282 = t312 * qJD(4) + t296 * t317;
t280 = t311 * qJD(1) + t286 * qJD(2);
t279 = t286 * qJD(1) + t311 * qJD(2);
t278 = -qJD(1) * t321 - t306 * t324 + t316 * t330;
t273 = t296 * t318 - qJD(4) * t334 + (qJD(4) * t332 - t279) * t297;
t272 = t313 * qJD(4) - t279 * t296 - t297 * t318;
t1 = [-t315 * qJD(5) - t281 * t295 - t285 * qJD(3) - t336 * t280 + t337 * t341 + t335 * t340 + (-t306 * pkin(1) - t320 * t332) * qJD(1), t288 * qJD(3) + t338 * t278 - t336 * t279 + t311 * t342, -t278, t313 * qJD(5) - t337 * t272 + t335 * t273, t272, 0; -(t297 * t332 - t334) * qJD(5) - t279 * t295 + t311 * qJD(3) - t336 * t278 + t337 * t273 + t335 * t272 + (-t304 * pkin(1) + t320 * t331) * qJD(1), t286 * qJD(3) - t280 * t338 + t336 * t281 + t342 * t285, t280, -t314 * qJD(5) - t335 * t341 + t337 * t340, -t340, 0; 0 (qJD(3) * t303 - t342 * t305 + (-t303 * t338 + t336 * t305) * qJD(2)) * t300, t300 * t325, t312 * qJD(5) + t335 * (t297 * t317 + (-t296 * t333 + t297 * t301) * qJD(4)) - t337 * t282, t282, 0;];
JaD_transl  = t1;
