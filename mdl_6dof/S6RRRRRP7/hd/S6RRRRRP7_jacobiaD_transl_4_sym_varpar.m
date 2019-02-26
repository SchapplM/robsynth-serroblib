% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP7_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP7_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:04
% EndTime: 2019-02-26 22:43:04
% DurationCPUTime: 0.30s
% Computational Cost: add. (332->69), mult. (663->114), div. (0->0), fcn. (618->10), ass. (0->55)
t279 = cos(pkin(6));
t281 = sin(qJ(2));
t282 = sin(qJ(1));
t313 = t282 * t281;
t301 = t279 * t313;
t284 = cos(qJ(2));
t285 = cos(qJ(1));
t310 = t285 * t284;
t265 = -qJD(1) * t301 - qJD(2) * t313 + (qJD(2) * t279 + qJD(1)) * t310;
t276 = qJD(3) + qJD(4);
t278 = sin(pkin(6));
t314 = t278 * t285;
t321 = t276 * t314 - t265;
t320 = r_i_i_C(3) + pkin(10) + pkin(9);
t319 = pkin(3) * qJD(3);
t311 = t285 * t281;
t312 = t282 * t284;
t267 = t279 * t311 + t312;
t318 = t267 * t276;
t277 = qJ(3) + qJ(4);
t274 = sin(t277);
t317 = t274 * t276;
t316 = t278 * t282;
t283 = cos(qJ(3));
t315 = t278 * t283;
t275 = cos(t277);
t291 = t301 - t310;
t305 = qJD(1) * t285;
t299 = t278 * t305;
t289 = t276 * t291 + t299;
t292 = t279 * t312 + t311;
t263 = t267 * qJD(1) + t292 * qJD(2);
t295 = t276 * t316 - t263;
t258 = -t295 * t274 + t289 * t275;
t259 = t289 * t274 + t295 * t275;
t309 = t258 * r_i_i_C(1) - t259 * r_i_i_C(2);
t306 = qJD(1) * t282;
t300 = t278 * t306;
t290 = t300 - t318;
t297 = t321 * t275;
t308 = (t321 * t274 + t290 * t275) * r_i_i_C(1) + (-t290 * t274 + t297) * r_i_i_C(2);
t298 = qJD(2) * t278 * t284;
t288 = -t276 * t279 - t298;
t303 = t276 * t278 * t281;
t307 = (t288 * t274 - t275 * t303) * r_i_i_C(1) + (t274 * t303 + t288 * t275) * r_i_i_C(2);
t280 = sin(qJ(3));
t304 = t280 * t319;
t296 = -r_i_i_C(1) * t274 - r_i_i_C(2) * t275;
t273 = t283 * pkin(3) + pkin(2);
t294 = t275 * r_i_i_C(1) - t274 * r_i_i_C(2) + t273;
t293 = t279 * t310 - t313;
t287 = t296 * t276 - t304;
t264 = t292 * qJD(1) + t267 * qJD(2);
t262 = -t293 * qJD(1) + t291 * qJD(2);
t1 = [(t267 * t317 + t297) * r_i_i_C(1) + (t265 * t274 + t275 * t318) * r_i_i_C(2) - t265 * t273 + t267 * t304 - pkin(1) * t305 - t320 * t264 + ((-r_i_i_C(2) * t317 + t283 * t319) * t285 + (-pkin(3) * t280 - pkin(8) + t296) * t306) * t278, t294 * t262 - t320 * t263 - t287 * t292 (t283 * t299 + t263 * t280 + (-t280 * t316 + t283 * t291) * qJD(3)) * pkin(3) + t309, t309, 0, 0; t259 * r_i_i_C(1) + t258 * r_i_i_C(2) - t263 * t273 - t320 * t262 + (-pkin(1) * t282 + pkin(8) * t314) * qJD(1) + (t280 * t299 + (t280 * t291 + t282 * t315) * qJD(3)) * pkin(3), -t294 * t264 + t320 * t265 + t287 * t293 (t283 * t300 - t265 * t280 + (-t267 * t283 + t280 * t314) * qJD(3)) * pkin(3) + t308, t308, 0, 0; 0 (t287 * t284 + (-t294 * t281 + t320 * t284) * qJD(2)) * t278 (-t280 * t298 + (-t279 * t280 - t281 * t315) * qJD(3)) * pkin(3) + t307, t307, 0, 0;];
JaD_transl  = t1;
