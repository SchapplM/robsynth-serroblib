% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR3
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiaD_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:24
% EndTime: 2019-02-26 20:05:25
% DurationCPUTime: 0.30s
% Computational Cost: add. (198->86), mult. (654->154), div. (0->0), fcn. (660->12), ass. (0->49)
t275 = sin(pkin(12));
t277 = sin(pkin(6));
t304 = t275 * t277;
t276 = sin(pkin(7));
t303 = t276 * t277;
t279 = cos(pkin(12));
t302 = t277 * t279;
t280 = cos(pkin(7));
t284 = cos(qJ(3));
t301 = t280 * t284;
t281 = cos(pkin(6));
t283 = sin(qJ(2));
t300 = t281 * t283;
t285 = cos(qJ(2));
t299 = t281 * t285;
t282 = sin(qJ(3));
t298 = t282 * t285;
t297 = t283 * t284;
t296 = qJD(2) * t283;
t295 = qJD(3) * t282;
t294 = qJD(3) * t284;
t293 = pkin(3) * t295;
t292 = t279 * t299;
t291 = t280 * t294;
t274 = sin(pkin(13));
t278 = cos(pkin(13));
t290 = t284 * t274 + t282 * t278;
t265 = t282 * t274 - t284 * t278;
t259 = t275 * t285 + t279 * t300;
t289 = t275 * t299 + t279 * t283;
t288 = t275 * t300 - t279 * t285;
t287 = qJD(3) * t290;
t262 = t265 * qJD(3);
t251 = t265 * t280;
t252 = t290 * t280;
t286 = t280 * t282 * pkin(3) + t252 * r_i_i_C(1) - t251 * r_i_i_C(2) + (-r_i_i_C(3) - pkin(9) - qJ(4)) * t276;
t273 = t284 * pkin(3) + pkin(2);
t264 = pkin(3) * t291 - t276 * qJD(4);
t263 = -t274 * t294 - t278 * t295;
t258 = -t275 * t283 + t292;
t256 = t288 * qJD(2);
t255 = t289 * qJD(2);
t254 = t259 * qJD(2);
t253 = -qJD(2) * t292 + t275 * t296;
t250 = t280 * t287;
t249 = t280 * t274 * t295 - t278 * t291;
t248 = t276 * t287;
t247 = t276 * t262;
t1 = [0 (-t249 * t288 - t256 * t265 - t263 * t289) * r_i_i_C(1) + (-t250 * t288 - t256 * t290 - t262 * t289) * r_i_i_C(2) + t256 * t273 + t289 * t293 + t288 * t264 + t286 * t255 (-t248 * t304 + t250 * t289 - t256 * t251 + t255 * t290 - t262 * t288) * r_i_i_C(1) + (t247 * t304 - t249 * t289 - t256 * t252 - t255 * t265 + t263 * t288) * r_i_i_C(2) + (t256 * t301 + t255 * t282 + (t288 * t284 + (-t275 * t303 + t280 * t289) * t282) * qJD(3)) * pkin(3), -t256 * t276, 0, 0; 0 (t259 * t249 + t254 * t265 + t258 * t263) * r_i_i_C(1) + (t259 * t250 + t254 * t290 + t258 * t262) * r_i_i_C(2) - t254 * t273 - t258 * t293 - t259 * t264 + t286 * t253 (t248 * t302 - t258 * t250 + t254 * t251 + t253 * t290 + t259 * t262) * r_i_i_C(1) + (-t247 * t302 + t258 * t249 + t254 * t252 - t253 * t265 - t259 * t263) * r_i_i_C(2) + (-t254 * t301 + t253 * t282 + (-t259 * t284 + (-t258 * t280 + t276 * t302) * t282) * qJD(3)) * pkin(3), t254 * t276, 0, 0; 0 ((t249 * t283 + t263 * t285) * r_i_i_C(1) + (t250 * t283 + t262 * t285) * r_i_i_C(2) - t285 * t293 - t283 * t264 + ((t265 * r_i_i_C(1) + r_i_i_C(2) * t290 - t273) * t283 - t286 * t285) * qJD(2)) * t277 (-t248 * r_i_i_C(1) + t247 * r_i_i_C(2) - t276 * t293) * t281 + ((-t250 * t285 + t262 * t283) * r_i_i_C(1) + (t249 * t285 - t263 * t283) * r_i_i_C(2) + ((t251 * t283 - t285 * t290) * r_i_i_C(1) + (t252 * t283 + t265 * t285) * r_i_i_C(2)) * qJD(2) + ((-t280 * t298 - t297) * qJD(3) + (-t280 * t297 - t298) * qJD(2)) * pkin(3)) * t277, t296 * t303, 0, 0;];
JaD_transl  = t1;
