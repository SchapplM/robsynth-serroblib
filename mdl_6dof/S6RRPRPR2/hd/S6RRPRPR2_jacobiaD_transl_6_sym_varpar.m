% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:17
% EndTime: 2019-02-26 21:38:17
% DurationCPUTime: 0.29s
% Computational Cost: add. (485->65), mult. (479->88), div. (0->0), fcn. (361->10), ass. (0->54)
t269 = qJ(2) + pkin(10);
t266 = qJ(4) + t269;
t263 = cos(t266);
t270 = sin(qJ(6));
t273 = cos(qJ(6));
t318 = r_i_i_C(1) * t270 + r_i_i_C(2) * t273 + qJ(5);
t322 = t263 * t318;
t258 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t269);
t247 = t258 * qJD(2);
t268 = qJD(2) + qJD(4);
t311 = t263 * t268;
t257 = qJ(5) * t311;
t262 = sin(t266);
t302 = pkin(4) + pkin(9) + r_i_i_C(3);
t292 = t302 * t268;
t321 = (qJD(5) - t292) * t262 + (pkin(5) + pkin(8) + qJ(3) + pkin(7)) * qJD(1) + t247 + t257;
t303 = qJD(6) * t273;
t295 = t263 * t303;
t320 = r_i_i_C(1) * t295 + qJD(5) * t263;
t290 = qJD(6) * t262 + qJD(1);
t317 = t273 * t290;
t316 = t290 * t270;
t310 = t268 * t270;
t309 = t268 * t273;
t275 = cos(qJ(1));
t308 = t268 * t275;
t272 = sin(qJ(1));
t307 = qJD(1) * t272;
t306 = qJD(1) * t275;
t304 = qJD(6) * t270;
t301 = t263 * t309;
t300 = t272 * t311;
t299 = t263 * t308;
t297 = t262 * t307;
t296 = t263 * t304;
t294 = t302 * t262;
t293 = t302 * t263;
t289 = -qJD(1) * t262 - qJD(6);
t288 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t269);
t286 = t272 * t320 + t306 * t322;
t285 = t275 * t320 + t302 * t297;
t284 = t289 * t275;
t282 = t318 * t262;
t281 = -r_i_i_C(2) * t304 - t292;
t280 = t289 * t272 + t299;
t279 = qJD(3) + (-qJ(5) * t262 - pkin(1) + t288 - t293) * qJD(1);
t278 = t263 * r_i_i_C(1) * t310 + r_i_i_C(2) * t301 + t257 + (r_i_i_C(1) * t303 + qJD(5) + t281) * t262;
t277 = -r_i_i_C(2) * t296 + (-t293 - t282) * t268;
t276 = t288 * qJD(2) + t277;
t242 = t280 * t270 + t275 * t317;
t241 = t280 * t273 - t275 * t316;
t240 = -t272 * t317 + (t284 - t300) * t270;
t239 = t273 * t284 + (-t301 + t316) * t272;
t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) - t321 * t272 + t279 * t275 (-t258 - t322) * t307 + t276 * t275 + t285, t306, -t282 * t308 + (t281 * t275 - t307 * t318) * t263 + t285, -t297 + t299, t241 * r_i_i_C(1) - t242 * r_i_i_C(2); t242 * r_i_i_C(1) + t241 * r_i_i_C(2) + t279 * t272 + t321 * t275 (t258 - t294) * t306 + t276 * t272 + t286, t307, t277 * t272 - t294 * t306 + t286, t262 * t306 + t300, -t239 * r_i_i_C(1) + t240 * r_i_i_C(2); 0, t247 + t278, 0, t278, t268 * t262 (-t262 * t310 + t295) * r_i_i_C(2) + (t262 * t309 + t296) * r_i_i_C(1);];
JaD_transl  = t1;
