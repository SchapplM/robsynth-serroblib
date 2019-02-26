% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPP3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:03
% EndTime: 2019-02-26 20:10:03
% DurationCPUTime: 0.38s
% Computational Cost: add. (224->70), mult. (725->134), div. (0->0), fcn. (724->10), ass. (0->50)
t290 = sin(qJ(3));
t293 = cos(qJ(3));
t289 = sin(qJ(4));
t292 = cos(qJ(4));
t304 = r_i_i_C(1) * t292 - r_i_i_C(2) * t289;
t302 = pkin(3) + t304;
t323 = pkin(9) + r_i_i_C(3);
t325 = (t302 * t290 - t323 * t293) * qJD(3);
t303 = r_i_i_C(1) * t289 + r_i_i_C(2) * t292;
t320 = cos(pkin(6));
t287 = sin(pkin(6));
t319 = t287 * t290;
t318 = t287 * t293;
t294 = cos(qJ(2));
t317 = t287 * t294;
t291 = sin(qJ(2));
t316 = qJD(2) * t291;
t315 = qJD(2) * t294;
t314 = qJD(4) * t289;
t313 = qJD(4) * t292;
t312 = qJD(4) * t293;
t311 = t287 * t316;
t310 = t287 * t315;
t309 = t291 * t320;
t308 = t294 * t320;
t286 = sin(pkin(10));
t306 = t286 * t309;
t288 = cos(pkin(10));
t305 = t288 * t308;
t278 = t286 * t294 + t288 * t309;
t301 = -t278 * t290 - t288 * t318;
t300 = -t278 * t293 + t288 * t319;
t280 = t288 * t294 - t306;
t299 = -t280 * t290 + t286 * t318;
t270 = t280 * t293 + t286 * t319;
t298 = qJD(4) * t303;
t279 = t286 * t308 + t288 * t291;
t297 = -t291 * t319 + t320 * t293;
t282 = t320 * t290 + t291 * t318;
t296 = -t323 * t290 - t302 * t293 - pkin(2);
t295 = t303 * t312 + t325;
t277 = t286 * t291 - t305;
t276 = -qJD(2) * t306 + t288 * t315;
t275 = t279 * qJD(2);
t274 = t278 * qJD(2);
t273 = -qJD(2) * t305 + t286 * t316;
t272 = t297 * qJD(3) + t293 * t310;
t266 = t299 * qJD(3) - t275 * t293;
t264 = t301 * qJD(3) - t273 * t293;
t1 = [0 (-t275 * t289 + t280 * t313) * r_i_i_C(1) + (-t275 * t292 - t280 * t314) * r_i_i_C(2) - t275 * pkin(8) + t296 * t276 + t295 * t279, t323 * t266 - t299 * t298 + t302 * (-t270 * qJD(3) + t275 * t290) (-t266 * t289 + t276 * t292) * r_i_i_C(1) + (-t266 * t292 - t276 * t289) * r_i_i_C(2) + ((-t270 * t292 - t279 * t289) * r_i_i_C(1) + (t270 * t289 - t279 * t292) * r_i_i_C(2)) * qJD(4), 0, 0; 0 (-t273 * t289 + t278 * t313) * r_i_i_C(1) + (-t273 * t292 - t278 * t314) * r_i_i_C(2) - t273 * pkin(8) + t296 * t274 + t295 * t277, t323 * t264 - t301 * t298 + t302 * (t300 * qJD(3) + t273 * t290) (-t264 * t289 + t274 * t292) * r_i_i_C(1) + (-t264 * t292 - t274 * t289) * r_i_i_C(2) + ((-t277 * t289 + t292 * t300) * r_i_i_C(1) + (-t277 * t292 - t289 * t300) * r_i_i_C(2)) * qJD(4), 0, 0; 0 ((t296 * qJD(2) + t304 * qJD(4)) * t291 + (qJD(2) * pkin(8) - t325 + t303 * (qJD(2) - t312)) * t294) * t287, t323 * t272 - t297 * t298 + t302 * (-t282 * qJD(3) - t290 * t310) (-t272 * t289 + t292 * t311) * r_i_i_C(1) + (-t272 * t292 - t289 * t311) * r_i_i_C(2) + ((-t282 * t292 + t289 * t317) * r_i_i_C(1) + (t282 * t289 + t292 * t317) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
