% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:43
% EndTime: 2019-02-26 20:00:43
% DurationCPUTime: 0.19s
% Computational Cost: add. (203->43), mult. (659->75), div. (0->0), fcn. (644->10), ass. (0->41)
t262 = sin(pkin(10));
t265 = cos(pkin(10));
t270 = cos(qJ(2));
t266 = cos(pkin(6));
t268 = sin(qJ(2));
t288 = t266 * t268;
t256 = t262 * t270 + t265 * t288;
t269 = cos(qJ(3));
t263 = sin(pkin(6));
t267 = sin(qJ(3));
t290 = t263 * t267;
t293 = -t256 * t269 + t265 * t290;
t261 = sin(pkin(11));
t264 = cos(pkin(11));
t281 = t261 * r_i_i_C(1) + t264 * r_i_i_C(2) + qJ(4);
t285 = pkin(3) + r_i_i_C(3) + qJ(5);
t292 = t281 * t267 + t285 * t269 + pkin(2);
t289 = t263 * t269;
t287 = t266 * t270;
t286 = qJD(2) * t268;
t283 = t265 * t287;
t282 = qJD(2) * t263 * t270;
t280 = t264 * r_i_i_C(1) - t261 * r_i_i_C(2) + pkin(4) + pkin(8);
t279 = t256 * t267 + t265 * t289;
t275 = t262 * t288 - t265 * t270;
t278 = t262 * t289 + t267 * t275;
t277 = t262 * t290 - t269 * t275;
t276 = t262 * t287 + t265 * t268;
t274 = t266 * t267 + t268 * t289;
t273 = -t266 * t269 + t268 * t290;
t272 = qJD(2) * t292;
t271 = t267 * qJD(4) + t269 * qJD(5) + (-t285 * t267 + t281 * t269) * qJD(3);
t253 = t276 * qJD(2);
t251 = -qJD(2) * t283 + t262 * t286;
t250 = -t273 * qJD(3) + t269 * t282;
t249 = t274 * qJD(3) + t267 * t282;
t248 = t278 * qJD(3) - t253 * t269;
t247 = t277 * qJD(3) - t253 * t267;
t246 = -t279 * qJD(3) - t251 * t269;
t245 = -t293 * qJD(3) - t251 * t267;
t1 = [0, -t280 * t253 - t271 * t276 + t275 * t272, t277 * qJD(4) + t278 * qJD(5) - t285 * t247 + t281 * t248, t247, t248, 0; 0, -t280 * t251 - t256 * t272 + t271 * (-t262 * t268 + t283) -t293 * qJD(4) - t279 * qJD(5) - t285 * t245 + t281 * t246, t245, t246, 0; 0 (-t292 * t286 + (t280 * qJD(2) + t271) * t270) * t263, t274 * qJD(4) - t273 * qJD(5) - t285 * t249 + t281 * t250, t249, t250, 0;];
JaD_transl  = t1;
