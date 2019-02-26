% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:01
% EndTime: 2019-02-26 20:00:01
% DurationCPUTime: 0.32s
% Computational Cost: add. (241->65), mult. (793->118), div. (0->0), fcn. (787->10), ass. (0->53)
t287 = sin(pkin(10));
t290 = cos(pkin(10));
t295 = cos(qJ(2));
t291 = cos(pkin(6));
t293 = sin(qJ(2));
t319 = t291 * t293;
t280 = t287 * t295 + t290 * t319;
t294 = cos(qJ(3));
t288 = sin(pkin(6));
t292 = sin(qJ(3));
t321 = t288 * t292;
t326 = -t280 * t294 + t290 * t321;
t325 = pkin(4) + r_i_i_C(1);
t324 = r_i_i_C(2) + qJ(4);
t323 = r_i_i_C(3) + qJ(5);
t320 = t288 * t294;
t318 = t291 * t295;
t317 = t293 * t294;
t316 = qJD(2) * t293;
t315 = qJD(2) * t295;
t314 = qJD(3) * t292;
t286 = sin(pkin(11));
t313 = qJD(5) * t286;
t289 = cos(pkin(11));
t312 = t289 * qJD(5);
t310 = t287 * t316;
t309 = t288 * t315;
t308 = t290 * t315;
t307 = t295 * t314;
t306 = -t280 * t292 - t290 * t320;
t282 = -t287 * t319 + t290 * t295;
t305 = -t282 * t292 + t287 * t320;
t304 = t282 * t294 + t287 * t321;
t303 = t287 * t318 + t290 * t293;
t302 = t288 * t317 + t291 * t292;
t301 = t291 * t294 - t293 * t321;
t276 = t280 * qJD(2);
t279 = -t287 * t293 + t290 * t318;
t300 = -t276 * t294 - t279 * t314;
t278 = -t291 * t310 + t308;
t299 = -t278 * t294 + t303 * t314;
t298 = -pkin(3) * t294 - t292 * t324 - pkin(2);
t297 = -t286 * t323 - t289 * t325 - pkin(3);
t296 = t294 * t313 + t292 * qJD(4) + (-pkin(3) * t292 + t294 * t324) * qJD(3);
t277 = t303 * qJD(2);
t275 = -t291 * t308 + t310;
t274 = qJD(3) * t301 + t294 * t309;
t273 = qJD(3) * t302 + t292 * t309;
t270 = qJD(3) * t305 - t277 * t294;
t269 = qJD(3) * t304 - t277 * t292;
t268 = qJD(3) * t306 - t275 * t294;
t267 = -qJD(3) * t326 - t275 * t292;
t1 = [0, -t282 * t312 - t277 * pkin(8) + t325 * (-t277 * t286 + t289 * t299) + t323 * (t277 * t289 + t286 * t299) - t296 * t303 + t298 * t278, qJD(4) * t304 + t269 * t297 + t270 * t324 + t305 * t313, t269, t270 * t286 - t278 * t289, 0; 0, -t280 * t312 - t275 * pkin(8) + t325 * (-t275 * t286 + t289 * t300) + t323 * (t275 * t289 + t286 * t300) + t296 * t279 + t298 * t276, -qJD(4) * t326 + t267 * t297 + t268 * t324 + t306 * t313, t267, t268 * t286 - t276 * t289, 0; 0 (t325 * (-t289 * t307 + (t286 * t295 - t289 * t317) * qJD(2)) - t323 * (t286 * t307 + (t286 * t317 + t289 * t295) * qJD(2)) - t293 * t312 + t296 * t295 + (t295 * pkin(8) + t298 * t293) * qJD(2)) * t288, qJD(4) * t302 + t273 * t297 + t274 * t324 + t301 * t313, t273, -t288 * t289 * t316 + t274 * t286, 0;];
JaD_transl  = t1;
