% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:35
% EndTime: 2019-02-26 19:51:35
% DurationCPUTime: 0.37s
% Computational Cost: add. (360->76), mult. (750->138), div. (0->0), fcn. (751->11), ass. (0->53)
t303 = pkin(11) + qJ(4);
t301 = sin(t303);
t302 = cos(t303);
t308 = sin(qJ(5));
t310 = cos(qJ(5));
t321 = r_i_i_C(1) * t310 - r_i_i_C(2) * t308;
t319 = pkin(4) + t321;
t341 = pkin(9) + r_i_i_C(3);
t343 = (t301 * t319 - t302 * t341) * qJD(4);
t320 = r_i_i_C(1) * t308 + r_i_i_C(2) * t310;
t338 = cos(pkin(6));
t304 = sin(pkin(10));
t305 = sin(pkin(6));
t337 = t304 * t305;
t306 = cos(pkin(10));
t336 = t305 * t306;
t309 = sin(qJ(2));
t335 = t305 * t309;
t311 = cos(qJ(2));
t334 = t305 * t311;
t333 = qJD(2) * t309;
t332 = qJD(2) * t311;
t331 = qJD(5) * t302;
t330 = qJD(5) * t308;
t329 = qJD(5) * t310;
t328 = t305 * t332;
t327 = t305 * t333;
t326 = t309 * t338;
t325 = t311 * t338;
t323 = t304 * t326;
t322 = t306 * t325;
t294 = t304 * t311 + t306 * t326;
t318 = -t294 * t301 - t302 * t336;
t317 = -t294 * t302 + t301 * t336;
t296 = t306 * t311 - t323;
t316 = -t296 * t301 + t302 * t337;
t286 = t296 * t302 + t301 * t337;
t315 = qJD(5) * t320;
t295 = t304 * t325 + t306 * t309;
t314 = -t301 * t335 + t302 * t338;
t288 = t301 * t338 + t302 * t335;
t313 = -t341 * t301 - t319 * t302 - cos(pkin(11)) * pkin(3) - pkin(2);
t312 = t320 * t331 + t343;
t307 = -pkin(8) - qJ(3);
t293 = t304 * t309 - t322;
t292 = -qJD(2) * t323 + t306 * t332;
t291 = t295 * qJD(2);
t290 = t294 * qJD(2);
t289 = -qJD(2) * t322 + t304 * t333;
t282 = qJD(4) * t314 + t302 * t328;
t280 = qJD(4) * t316 - t291 * t302;
t278 = qJD(4) * t318 - t289 * t302;
t1 = [0 (-t291 * t308 + t296 * t329) * r_i_i_C(1) + (-t291 * t310 - t296 * t330) * r_i_i_C(2) + t291 * t307 + t296 * qJD(3) + t313 * t292 + t312 * t295, t292, t341 * t280 - t316 * t315 + t319 * (-qJD(4) * t286 + t291 * t301) (-t280 * t308 + t292 * t310) * r_i_i_C(1) + (-t280 * t310 - t292 * t308) * r_i_i_C(2) + ((-t286 * t310 - t295 * t308) * r_i_i_C(1) + (t286 * t308 - t295 * t310) * r_i_i_C(2)) * qJD(5), 0; 0 (-t289 * t308 + t294 * t329) * r_i_i_C(1) + (-t289 * t310 - t294 * t330) * r_i_i_C(2) + t289 * t307 + t294 * qJD(3) + t313 * t290 + t312 * t293, t290, t341 * t278 - t318 * t315 + t319 * (qJD(4) * t317 + t289 * t301) (-t278 * t308 + t290 * t310) * r_i_i_C(1) + (-t278 * t310 - t290 * t308) * r_i_i_C(2) + ((-t293 * t308 + t310 * t317) * r_i_i_C(1) + (-t293 * t310 - t308 * t317) * r_i_i_C(2)) * qJD(5), 0; 0 ((qJD(2) * t313 + qJD(5) * t321 + qJD(3)) * t309 + (-qJD(2) * t307 - t343 + t320 * (qJD(2) - t331)) * t311) * t305, t327, t341 * t282 - t314 * t315 + t319 * (-qJD(4) * t288 - t301 * t328) (-t282 * t308 + t310 * t327) * r_i_i_C(1) + (-t282 * t310 - t308 * t327) * r_i_i_C(2) + ((-t288 * t310 + t308 * t334) * r_i_i_C(1) + (t288 * t308 + t310 * t334) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
