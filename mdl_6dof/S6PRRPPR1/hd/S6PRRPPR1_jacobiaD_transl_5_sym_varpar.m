% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:09
% EndTime: 2019-02-26 19:58:10
% DurationCPUTime: 0.21s
% Computational Cost: add. (276->55), mult. (590->97), div. (0->0), fcn. (568->12), ass. (0->43)
t277 = sin(pkin(10));
t280 = cos(pkin(10));
t286 = cos(qJ(2));
t281 = cos(pkin(6));
t284 = sin(qJ(2));
t301 = t281 * t284;
t266 = t277 * t286 + t280 * t301;
t275 = qJ(3) + pkin(11);
t273 = sin(t275);
t274 = cos(t275);
t278 = sin(pkin(6));
t304 = t278 * t280;
t308 = -t266 * t274 + t273 * t304;
t307 = r_i_i_C(3) + qJ(5);
t305 = t277 * t278;
t283 = sin(qJ(3));
t303 = t278 * t283;
t302 = t278 * t284;
t300 = t281 * t286;
t299 = qJD(2) * t284;
t298 = qJD(2) * t286;
t296 = t277 * t299;
t295 = t278 * t298;
t294 = t280 * t298;
t276 = sin(pkin(12));
t279 = cos(pkin(12));
t293 = -r_i_i_C(1) * t279 + r_i_i_C(2) * t276 - pkin(4);
t292 = r_i_i_C(1) * t276 + r_i_i_C(2) * t279 + pkin(8) + qJ(4);
t268 = -t277 * t301 + t280 * t286;
t291 = t268 * t274 + t273 * t305;
t290 = t273 * t281 + t274 * t302;
t289 = t277 * t300 + t280 * t284;
t285 = cos(qJ(3));
t288 = -pkin(3) * t285 - t307 * t273 + t293 * t274 - pkin(2);
t287 = t273 * qJD(5) + (-pkin(3) * t283 + t293 * t273 + t307 * t274) * qJD(3);
t264 = -t281 * t296 + t294;
t263 = t289 * qJD(2);
t262 = t266 * qJD(2);
t261 = -t281 * t294 + t296;
t259 = t290 * qJD(3) + t273 * t295;
t257 = t291 * qJD(3) - t263 * t273;
t255 = -t308 * qJD(3) - t261 * t273;
t1 = [0, t268 * qJD(4) - t292 * t263 + t288 * t264 - t287 * t289, t291 * qJD(5) + t307 * (-t263 * t274 + (-t268 * t273 + t274 * t305) * qJD(3)) + t293 * t257 + (t263 * t283 + (-t268 * t285 - t277 * t303) * qJD(3)) * pkin(3), t264, t257, 0; 0, t266 * qJD(4) - t292 * t261 + t288 * t262 + t287 * (-t277 * t284 + t280 * t300) -t308 * qJD(5) + t307 * (-t261 * t274 + (-t266 * t273 - t274 * t304) * qJD(3)) + t293 * t255 + (t261 * t283 + (-t266 * t285 + t280 * t303) * qJD(3)) * pkin(3), t262, t255, 0; 0 ((t288 * qJD(2) + qJD(4)) * t284 + (t292 * qJD(2) + t287) * t286) * t278, t290 * qJD(5) + t307 * (t274 * t295 + (-t273 * t302 + t274 * t281) * qJD(3)) + t293 * t259 + (-t283 * t295 + (-t281 * t283 - t285 * t302) * qJD(3)) * pkin(3), t278 * t299, t259, 0;];
JaD_transl  = t1;
