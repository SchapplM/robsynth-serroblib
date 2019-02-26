% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:38
% EndTime: 2019-02-26 20:07:39
% DurationCPUTime: 0.35s
% Computational Cost: add. (275->63), mult. (880->119), div. (0->0), fcn. (880->10), ass. (0->47)
t301 = sin(qJ(3));
t304 = cos(qJ(3));
t300 = sin(qJ(5));
t303 = cos(qJ(5));
t319 = t303 * r_i_i_C(1) - t300 * r_i_i_C(2);
t308 = t319 * qJD(5) + qJD(4);
t318 = -t300 * r_i_i_C(1) - t303 * r_i_i_C(2);
t316 = qJ(4) - t318;
t327 = pkin(3) + pkin(9) + r_i_i_C(3);
t306 = (t327 * t301 - t316 * t304) * qJD(3) - t308 * t301;
t298 = sin(pkin(11));
t302 = sin(qJ(2));
t305 = cos(qJ(2));
t333 = cos(pkin(11));
t334 = cos(pkin(6));
t317 = t334 * t333;
t289 = t298 * t305 + t302 * t317;
t299 = sin(pkin(6));
t323 = t299 * t333;
t336 = t289 * t304 - t301 * t323;
t324 = t298 * t334;
t291 = -t302 * t324 + t333 * t305;
t331 = t299 * t301;
t330 = t299 * t304;
t329 = t299 * t305;
t328 = qJD(2) * t302;
t326 = t299 * t328;
t325 = qJD(2) * t329;
t315 = t305 * t317;
t314 = pkin(4) + pkin(8) + t319;
t313 = -t291 * t301 + t298 * t330;
t312 = t291 * t304 + t298 * t331;
t311 = t318 * qJD(5);
t292 = t302 * t331 - t334 * t304;
t310 = t334 * t301 + t302 * t330;
t309 = -t289 * t301 - t304 * t323;
t290 = t333 * t302 + t305 * t324;
t307 = -t316 * t301 - t327 * t304 - pkin(2);
t288 = t298 * t302 - t315;
t287 = t291 * qJD(2);
t286 = t290 * qJD(2);
t285 = t289 * qJD(2);
t284 = -qJD(2) * t315 + t298 * t328;
t282 = t310 * qJD(3) + t301 * t325;
t276 = t312 * qJD(3) - t286 * t301;
t274 = t336 * qJD(3) - t284 * t301;
t1 = [0, -t314 * t286 + t307 * t287 + t306 * t290 + t291 * t311, t308 * t312 + t316 * (t313 * qJD(3) - t286 * t304) - t327 * t276, t276 (t276 * t303 - t287 * t300) * r_i_i_C(1) + (-t276 * t300 - t287 * t303) * r_i_i_C(2) + ((-t290 * t303 + t300 * t313) * r_i_i_C(1) + (t290 * t300 + t303 * t313) * r_i_i_C(2)) * qJD(5), 0; 0, -t314 * t284 + t307 * t285 + t306 * t288 + t289 * t311, t308 * t336 + t316 * (t309 * qJD(3) - t284 * t304) - t327 * t274, t274 (t274 * t303 - t285 * t300) * r_i_i_C(1) + (-t274 * t300 - t285 * t303) * r_i_i_C(2) + ((-t288 * t303 + t300 * t309) * r_i_i_C(1) + (t288 * t300 + t303 * t309) * r_i_i_C(2)) * qJD(5), 0; 0 ((t307 * qJD(2) + t311) * t302 + (t314 * qJD(2) - t306) * t305) * t299, t308 * t310 + t316 * (-t292 * qJD(3) + t304 * t325) - t327 * t282, t282 (t282 * t303 - t300 * t326) * r_i_i_C(1) + (-t282 * t300 - t303 * t326) * r_i_i_C(2) + ((-t292 * t300 + t303 * t329) * r_i_i_C(1) + (-t292 * t303 - t300 * t329) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
