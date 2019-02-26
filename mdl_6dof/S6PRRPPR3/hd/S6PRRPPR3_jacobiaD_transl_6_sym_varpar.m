% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:25
% EndTime: 2019-02-26 19:59:25
% DurationCPUTime: 0.31s
% Computational Cost: add. (329->66), mult. (1048->118), div. (0->0), fcn. (1040->10), ass. (0->46)
t295 = sin(pkin(10));
t297 = cos(pkin(10));
t304 = cos(qJ(2));
t298 = cos(pkin(6));
t301 = sin(qJ(2));
t325 = t298 * t301;
t287 = t295 * t304 + t297 * t325;
t303 = cos(qJ(3));
t296 = sin(pkin(6));
t300 = sin(qJ(3));
t328 = t296 * t300;
t330 = t287 * t303 - t297 * t328;
t299 = sin(qJ(6));
t302 = cos(qJ(6));
t317 = -t302 * r_i_i_C(1) + t299 * r_i_i_C(2);
t310 = pkin(5) + qJ(4) - t317;
t322 = pkin(3) + pkin(4) + pkin(9) + r_i_i_C(3);
t306 = t310 * t300 + t322 * t303 + pkin(2);
t327 = t296 * t303;
t326 = t296 * t304;
t324 = t298 * t304;
t323 = qJD(2) * t301;
t320 = t297 * t324;
t319 = qJD(2) * t326;
t318 = t296 * t323;
t316 = -t299 * r_i_i_C(1) - t302 * r_i_i_C(2);
t276 = t287 * t300 + t297 * t327;
t312 = t295 * t325 - t297 * t304;
t315 = t295 * t327 + t300 * t312;
t314 = t295 * t328 - t303 * t312;
t313 = t295 * t324 + t297 * t301;
t311 = t298 * t300 + t301 * t327;
t290 = -t298 * t303 + t301 * t328;
t309 = -pkin(8) + qJ(5) - t316;
t308 = t316 * qJD(6) + qJD(4);
t307 = t317 * qJD(6) - qJD(5);
t305 = t308 * t300 + (-t322 * t300 + t310 * t303) * qJD(3);
t286 = -t295 * t301 + t320;
t285 = t312 * qJD(2);
t284 = t313 * qJD(2);
t283 = t287 * qJD(2);
t282 = -qJD(2) * t320 + t295 * t323;
t280 = t311 * qJD(3) + t300 * t319;
t274 = t314 * qJD(3) - t284 * t300;
t272 = qJD(3) * t330 - t282 * t300;
t1 = [0, t309 * t284 + t306 * t285 - t305 * t313 - t307 * t312, t308 * t314 + t310 * (t315 * qJD(3) - t284 * t303) - t322 * t274, t274, t285 (-t274 * t299 + t285 * t302) * r_i_i_C(1) + (-t274 * t302 - t285 * t299) * r_i_i_C(2) + ((t299 * t313 + t302 * t315) * r_i_i_C(1) + (-t299 * t315 + t302 * t313) * r_i_i_C(2)) * qJD(6); 0, t309 * t282 - t306 * t283 + t305 * t286 + t307 * t287, t308 * t330 + t310 * (-t276 * qJD(3) - t282 * t303) - t322 * t272, t272, -t283 (-t272 * t299 - t283 * t302) * r_i_i_C(1) + (-t272 * t302 + t283 * t299) * r_i_i_C(2) + ((-t276 * t302 - t286 * t299) * r_i_i_C(1) + (t276 * t299 - t286 * t302) * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t306 + t307) * t301 + (-t309 * qJD(2) + t305) * t304) * t296, t308 * t311 + t310 * (-t290 * qJD(3) + t303 * t319) - t322 * t280, t280, -t318 (-t280 * t299 - t302 * t318) * r_i_i_C(1) + (-t280 * t302 + t299 * t318) * r_i_i_C(2) + ((-t290 * t302 - t299 * t326) * r_i_i_C(1) + (t290 * t299 - t302 * t326) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
