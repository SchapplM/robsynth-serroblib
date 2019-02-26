% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:13
% EndTime: 2019-02-26 19:49:14
% DurationCPUTime: 0.19s
% Computational Cost: add. (164->42), mult. (541->75), div. (0->0), fcn. (524->10), ass. (0->37)
t264 = cos(pkin(6));
t265 = sin(qJ(4));
t267 = cos(qJ(4));
t261 = sin(pkin(6));
t268 = cos(qJ(2));
t285 = t261 * t268;
t290 = -t264 * t267 + t265 * t285;
t259 = sin(pkin(11));
t262 = cos(pkin(11));
t275 = r_i_i_C(1) * t262 - r_i_i_C(2) * t259 + pkin(4);
t288 = r_i_i_C(3) + qJ(5);
t289 = t275 * t265 - t288 * t267 + qJ(3);
t287 = t261 * t265;
t286 = t261 * t267;
t266 = sin(qJ(2));
t284 = t264 * t266;
t282 = t264 * t268;
t281 = qJD(2) * t266;
t280 = qJD(2) * t268;
t260 = sin(pkin(10));
t278 = t260 * t281;
t263 = cos(pkin(10));
t277 = t263 * t280;
t276 = t261 * t281;
t274 = -r_i_i_C(1) * t259 - r_i_i_C(2) * t262 - pkin(2) - pkin(8);
t251 = t260 * t266 - t263 * t282;
t273 = -t251 * t265 + t263 * t286;
t253 = t260 * t282 + t263 * t266;
t272 = t253 * t265 + t260 * t286;
t271 = t260 * t268 + t263 * t284;
t269 = -t267 * qJD(5) + qJD(3) + (t288 * t265 + t275 * t267) * qJD(4);
t250 = -t264 * t278 + t277;
t248 = t271 * qJD(2);
t245 = -t290 * qJD(4) - t267 * t276;
t243 = t273 * qJD(4) + t248 * t267;
t241 = t272 * qJD(4) - t250 * t267;
t1 = [0, t274 * t250 - t289 * t253 * qJD(2) + t269 * (-t260 * t284 + t263 * t268) t250, t272 * qJD(5) - t288 * (-t250 * t265 + (-t253 * t267 + t260 * t287) * qJD(4)) - t275 * t241, t241, 0; 0, t274 * t248 - t289 * (-t264 * t277 + t278) + t269 * t271, t248, -t273 * qJD(5) + t288 * (t248 * t265 + (t251 * t267 + t263 * t287) * qJD(4)) + t275 * t243, -t243, 0; 0 (t289 * t280 + (t274 * qJD(2) + t269) * t266) * t261, t276, -t290 * qJD(5) - t288 * (-t265 * t276 + (t264 * t265 + t267 * t285) * qJD(4)) - t275 * t245, t245, 0;];
JaD_transl  = t1;
