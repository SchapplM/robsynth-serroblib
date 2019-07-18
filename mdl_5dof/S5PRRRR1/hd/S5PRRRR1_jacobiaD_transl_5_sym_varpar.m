% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:15
% EndTime: 2019-07-18 13:29:15
% DurationCPUTime: 0.19s
% Computational Cost: add. (192->49), mult. (310->87), div. (0->0), fcn. (243->8), ass. (0->47)
t252 = qJ(3) + qJ(4);
t249 = sin(t252);
t250 = cos(t252);
t251 = qJD(3) + qJD(4);
t253 = sin(qJ(5));
t256 = cos(qJ(5));
t282 = qJD(5) * t256;
t295 = t250 * t251 * t253 + t249 * t282;
t283 = qJD(5) * t253;
t274 = t249 * t283;
t294 = r_i_i_C(1) * t274 + t295 * r_i_i_C(2);
t254 = sin(qJ(3));
t293 = pkin(2) * t254;
t292 = r_i_i_C(1) * t256;
t291 = r_i_i_C(3) * t250;
t255 = sin(qJ(2));
t290 = t251 * t255;
t289 = t251 * t256;
t258 = cos(qJ(2));
t288 = t251 * t258;
t287 = t256 * t258;
t286 = qJD(2) * t255;
t285 = qJD(2) * t258;
t284 = qJD(3) * t254;
t281 = r_i_i_C(2) * t249 * t253;
t280 = t249 * t290;
t279 = t249 * t289;
t278 = t249 * t288;
t276 = t250 * t288;
t275 = t249 * t286;
t270 = qJD(2) * t281;
t268 = qJD(5) * t250 - qJD(2);
t267 = qJD(2) * t250 - qJD(5);
t266 = t294 * t258 + t275 * t292;
t265 = t294 * t255 + t258 * t270 + t285 * t291;
t264 = t268 * t253;
t263 = -t281 - t291;
t262 = -t249 * t285 - t250 * t290;
t261 = t267 * t255 + t278;
t260 = t250 * r_i_i_C(2) * t282 + t263 * t251 + (t250 * t283 + t279) * r_i_i_C(1);
t257 = cos(qJ(3));
t259 = -pkin(2) * qJD(3) * t257 + (-r_i_i_C(3) * t249 - t250 * t292) * t251;
t236 = -t267 * t287 + (t264 + t279) * t255;
t235 = t268 * t256 * t255 + (t267 * t258 - t280) * t253;
t234 = t261 * t256 + t258 * t264;
t233 = t261 * t253 - t268 * t287;
t1 = [0, t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t262 * r_i_i_C(3) + (t255 * t284 - t257 * t285) * pkin(2), t259 * t258 + (t263 + t293) * t286 + t266, -t276 * t292 - t255 * t270 + (-t250 * t286 - t278) * r_i_i_C(3) + t266, t233 * r_i_i_C(1) + t234 * r_i_i_C(2); 0, 0, pkin(2) * t284 + t260, t260, (t250 * t289 - t274) * r_i_i_C(2) + t295 * r_i_i_C(1); 0, -t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + (-t275 + t276) * r_i_i_C(3) + (-t257 * t286 - t258 * t284) * pkin(2), (-t249 * t292 - t293) * t285 + t259 * t255 + t265, -r_i_i_C(3) * t280 + t262 * t292 + t265, -t235 * r_i_i_C(1) + t236 * r_i_i_C(2);];
JaD_transl  = t1;
