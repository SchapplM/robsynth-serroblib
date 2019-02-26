% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:45
% EndTime: 2019-02-26 19:47:45
% DurationCPUTime: 0.21s
% Computational Cost: add. (244->48), mult. (780->94), div. (0->0), fcn. (810->10), ass. (0->44)
t293 = sin(pkin(11));
t296 = cos(pkin(11));
t302 = cos(qJ(2));
t314 = qJD(2) * t302;
t300 = sin(qJ(2));
t315 = qJD(2) * t300;
t323 = t293 * t315 - t296 * t314;
t322 = pkin(4) - r_i_i_C(2);
t321 = -pkin(8) - r_i_i_C(1);
t320 = r_i_i_C(3) + qJ(5);
t319 = pkin(2) * qJD(2);
t295 = sin(pkin(6));
t299 = sin(qJ(4));
t318 = t295 * t299;
t301 = cos(qJ(4));
t317 = t295 * t301;
t298 = cos(pkin(6));
t316 = t298 * t300;
t279 = t323 * t298;
t286 = -t293 * t314 - t296 * t315;
t294 = sin(pkin(10));
t297 = cos(pkin(10));
t311 = t297 * t279 - t294 * t286;
t271 = t294 * t279 + t297 * t286;
t309 = t302 * t293 + t300 * t296;
t282 = t309 * t295;
t310 = t282 * t301 + t298 * t299;
t308 = t300 * t293 - t302 * t296;
t284 = t309 * t298;
t274 = t297 * t284 - t294 * t308;
t307 = -t274 * t301 + t297 * t318;
t276 = -t294 * t284 - t297 * t308;
t306 = t276 * t301 + t294 * t318;
t305 = qJD(2) * t309;
t304 = t320 * t299 + t322 * t301 + pkin(3);
t303 = qJD(5) * t299 + (-t322 * t299 + t320 * t301) * qJD(4);
t285 = t308 * qJD(2);
t283 = t308 * t298;
t280 = t298 * t305;
t277 = t323 * t295;
t265 = t310 * qJD(4) - t277 * t299;
t263 = t306 * qJD(4) + t271 * t299;
t261 = -t307 * qJD(4) - t299 * t311;
t1 = [0, -t321 * t271 + (t294 * t316 - t297 * t302) * t319 + t303 * (t294 * t283 - t297 * t309) + t304 * (t294 * t280 + t297 * t285) 0, t306 * qJD(5) + t320 * (t271 * t301 + (-t276 * t299 + t294 * t317) * qJD(4)) - t322 * t263, t263, 0; 0, t321 * t311 + (-t294 * t302 - t297 * t316) * t319 + t303 * (-t297 * t283 - t294 * t309) + t304 * (-t297 * t280 + t294 * t285) 0, -t307 * qJD(5) + t320 * (-t311 * t301 + (-t274 * t299 - t297 * t317) * qJD(4)) - t322 * t261, t261, 0; 0, t321 * t277 + (-pkin(2) * t315 - t303 * t308 - t304 * t305) * t295, 0, t310 * qJD(5) + t320 * (-t277 * t301 + (-t282 * t299 + t298 * t301) * qJD(4)) - t322 * t265, t265, 0;];
JaD_transl  = t1;
