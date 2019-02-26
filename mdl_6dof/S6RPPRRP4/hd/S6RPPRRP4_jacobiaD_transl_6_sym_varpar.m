% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:06
% EndTime: 2019-02-26 20:32:06
% DurationCPUTime: 0.29s
% Computational Cost: add. (324->60), mult. (890->94), div. (0->0), fcn. (893->8), ass. (0->44)
t292 = sin(qJ(4));
t294 = cos(qJ(4));
t325 = pkin(8) + r_i_i_C(2);
t331 = t325 * t294;
t335 = (-pkin(4) * t292 + t331) * qJD(4);
t310 = t325 * t292;
t334 = t294 * pkin(4) + pkin(3) + t310;
t291 = sin(qJ(5));
t293 = cos(qJ(5));
t320 = r_i_i_C(3) + qJ(6);
t324 = r_i_i_C(1) + pkin(5);
t300 = t291 * t320 + t293 * t324;
t298 = pkin(4) + t300;
t332 = t292 * t298 - t331;
t318 = sin(pkin(9));
t319 = cos(pkin(9));
t322 = sin(qJ(1));
t323 = cos(qJ(1));
t286 = t318 * t323 - t319 * t322;
t330 = qJD(1) * t322;
t329 = qJD(1) * t323;
t299 = t291 * t324 - t293 * t320;
t328 = qJD(5) * t299 - qJD(6) * t291;
t327 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t285 = -t318 * t322 - t319 * t323;
t283 = t285 * qJD(1);
t315 = qJD(4) * t292;
t306 = t286 * t315;
t297 = -qJD(5) * t285 - t283 * t294 + t306;
t284 = t286 * qJD(1);
t313 = qJD(5) * t294;
t303 = t286 * t313 + t284;
t326 = t291 * t303 + t293 * t297;
t317 = t291 * t294;
t316 = t293 * t294;
t314 = qJD(4) * t294;
t307 = t285 * t315;
t302 = t285 * t291 + t286 * t316;
t301 = t285 * t316 - t286 * t291;
t295 = -t328 * t292 + (t294 * t298 + t310) * qJD(4);
t278 = (t285 * t313 + t283) * t291 + (qJD(5) * t286 + t284 * t294 + t307) * t293;
t277 = -qJD(5) * t301 - t283 * t293 + t284 * t317 + t291 * t307;
t273 = -qJD(5) * t302 - t283 * t317 - t284 * t293 + t291 * t306;
t1 = [-(t285 * t293 - t286 * t317) * qJD(6) - t284 * pkin(7) - t324 * t326 + t320 * (-t291 * t297 + t293 * t303) + t286 * t335 - qJ(2) * t330 + t334 * t283 + t327 * t323, t329, 0, -t284 * t332 + t285 * t295, -qJD(6) * t301 - t277 * t324 + t278 * t320, t277; -(t285 * t317 + t286 * t293) * qJD(6) + t283 * pkin(7) + t324 * t278 + t320 * t277 - t285 * t335 + qJ(2) * t329 + t334 * t284 + t327 * t322, t330, 0, t283 * t332 + t286 * t295, -qJD(6) * t302 - t273 * t324 + t320 * t326, t273; 0, 0, 0, qJD(4) * t332 + t294 * t328, t299 * t314 + (qJD(5) * t300 - qJD(6) * t293) * t292, -qJD(5) * t292 * t293 - t291 * t314;];
JaD_transl  = t1;
