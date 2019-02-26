% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:13
% EndTime: 2019-02-26 19:49:14
% DurationCPUTime: 0.32s
% Computational Cost: add. (348->75), mult. (856->131), div. (0->0), fcn. (860->12), ass. (0->54)
t306 = pkin(11) + qJ(6);
t304 = sin(t306);
t305 = cos(t306);
t326 = r_i_i_C(1) * t305 - r_i_i_C(2) * t304;
t321 = t326 * qJD(6);
t313 = sin(qJ(4));
t315 = cos(qJ(4));
t324 = cos(pkin(11)) * pkin(5) + pkin(4) + t326;
t341 = r_i_i_C(3) + pkin(9) + qJ(5);
t342 = t324 * t313 - t341 * t315 + qJ(3);
t309 = sin(pkin(6));
t340 = t309 * t313;
t314 = sin(qJ(2));
t339 = t309 * t314;
t338 = t309 * t315;
t316 = cos(qJ(2));
t337 = t309 * t316;
t311 = cos(pkin(6));
t336 = t311 * t314;
t335 = t311 * t316;
t334 = qJD(2) * t314;
t333 = qJD(2) * t316;
t332 = qJD(4) * t315;
t308 = sin(pkin(10));
t331 = t308 * t334;
t330 = t309 * t333;
t310 = cos(pkin(10));
t329 = t310 * t333;
t328 = qJD(4) * t340;
t327 = t309 * t334;
t325 = -r_i_i_C(1) * t304 - r_i_i_C(2) * t305;
t291 = t308 * t314 - t310 * t335;
t280 = -t291 * t313 + t310 * t338;
t323 = t291 * t315 + t310 * t340;
t293 = t308 * t335 + t310 * t314;
t278 = t293 * t313 + t308 * t338;
t292 = t308 * t316 + t310 * t336;
t322 = t311 * t313 + t315 * t337;
t320 = qJD(6) * t325;
t319 = -pkin(5) * sin(pkin(11)) - pkin(2) - pkin(8) + t325;
t317 = -t315 * qJD(5) + qJD(3) + t313 * t320 + (t341 * t313 + t324 * t315) * qJD(4);
t296 = t311 * t315 - t313 * t337;
t294 = -t308 * t336 + t310 * t316;
t290 = -t311 * t331 + t329;
t289 = t293 * qJD(2);
t288 = t292 * qJD(2);
t287 = -t311 * t329 + t331;
t282 = t311 * t332 - t315 * t327 - t316 * t328;
t281 = t322 * qJD(4) - t313 * t327;
t276 = t280 * qJD(4) + t288 * t315;
t275 = t323 * qJD(4) + t288 * t313;
t274 = t278 * qJD(4) - t290 * t315;
t273 = -t290 * t313 - t293 * t332 + t308 * t328;
t1 = [0, -t289 * t342 + t319 * t290 - t293 * t321 + t317 * t294, t290, t278 * qJD(5) - t341 * t273 + (t293 * t315 - t308 * t340) * t320 - t324 * t274, t274 (t273 * t304 - t289 * t305) * r_i_i_C(1) + (t273 * t305 + t289 * t304) * r_i_i_C(2) + ((-t278 * t305 - t294 * t304) * r_i_i_C(1) + (t278 * t304 - t294 * t305) * r_i_i_C(2)) * qJD(6); 0, -t287 * t342 + t319 * t288 - t291 * t321 + t317 * t292, t288, -t280 * qJD(5) + t341 * t275 + t324 * t276 + t323 * t320, -t276 (-t275 * t304 - t287 * t305) * r_i_i_C(1) + (-t275 * t305 + t287 * t304) * r_i_i_C(2) + ((t280 * t305 - t292 * t304) * r_i_i_C(1) + (-t280 * t304 - t292 * t305) * r_i_i_C(2)) * qJD(6); 0 ((t342 * qJD(2) + t321) * t316 + (t319 * qJD(2) + t317) * t314) * t309, t327, t296 * qJD(5) - t341 * t281 - t324 * t282 - t322 * t320, t282 (t281 * t304 + t305 * t330) * r_i_i_C(1) + (t281 * t305 - t304 * t330) * r_i_i_C(2) + ((-t296 * t305 - t304 * t339) * r_i_i_C(1) + (t296 * t304 - t305 * t339) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
