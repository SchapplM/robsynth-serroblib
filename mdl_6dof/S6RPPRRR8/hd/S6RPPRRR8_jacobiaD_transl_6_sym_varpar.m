% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:38:40
% EndTime: 2019-02-26 20:38:40
% DurationCPUTime: 0.25s
% Computational Cost: add. (372->58), mult. (438->90), div. (0->0), fcn. (346->9), ass. (0->52)
t244 = qJD(5) + qJD(6);
t248 = sin(qJ(5));
t286 = pkin(5) * t248;
t271 = qJD(5) * t286;
t245 = qJ(5) + qJ(6);
t241 = sin(t245);
t242 = cos(t245);
t294 = r_i_i_C(1) * t241 + r_i_i_C(2) * t242;
t254 = t294 * t244 + t271;
t281 = r_i_i_C(3) + pkin(9) + pkin(8);
t295 = -t281 * qJD(4) + t254;
t293 = -r_i_i_C(1) * t242 + r_i_i_C(2) * t241;
t243 = pkin(10) + qJ(4);
t238 = sin(t243);
t250 = cos(qJ(5));
t292 = t250 * (qJD(5) * t238 + qJD(1));
t291 = t281 * t238;
t237 = pkin(5) * t250 + pkin(4);
t239 = cos(t243);
t289 = t281 * t239 - pkin(3) * sin(pkin(10)) - t237 * t238 - qJ(2);
t249 = sin(qJ(1));
t277 = qJD(1) * t238;
t266 = t244 + t277;
t251 = cos(qJ(1));
t273 = qJD(4) * t251;
t269 = t239 * t273;
t288 = t266 * t249 - t269;
t274 = qJD(4) * t249;
t270 = t239 * t274;
t287 = t266 * t251 + t270;
t267 = -t238 * t244 - qJD(1);
t259 = t267 * t251;
t230 = t288 * t241 + t242 * t259;
t231 = t241 * t259 - t288 * t242;
t279 = -t230 * r_i_i_C(1) + t231 * r_i_i_C(2);
t260 = t267 * t249;
t232 = -t287 * t241 + t242 * t260;
t233 = t241 * t260 + t287 * t242;
t278 = t232 * r_i_i_C(1) - t233 * r_i_i_C(2);
t276 = qJD(1) * t249;
t275 = qJD(4) * t238;
t272 = qJD(5) * t250;
t268 = qJD(1) * t281;
t264 = -qJD(5) - t277;
t263 = -pkin(1) - pkin(7) - qJ(3) - t286;
t262 = pkin(5) * t272 + qJD(3);
t258 = t237 - t293;
t256 = qJD(1) * t258;
t255 = t293 * t239 * t244 + t294 * t275;
t253 = -t238 * t271 + qJD(2) + (t237 * t239 + t291) * qJD(4);
t240 = qJD(1) * t251;
t1 = [t231 * r_i_i_C(1) + t230 * r_i_i_C(2) - t262 * t249 + t253 * t251 + (t289 * t249 + t263 * t251) * qJD(1), t240, -t276 (t251 * t268 - t258 * t274) * t238 + (-t295 * t249 + t251 * t256) * t239 (-t249 * t292 + (t264 * t251 - t270) * t248) * pkin(5) + t278, t278; t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t262 * t251 + t253 * t249 + (t263 * t249 - t289 * t251) * qJD(1), t276, t240 (t249 * t268 + t258 * t273) * t238 + (t249 * t256 + t295 * t251) * t239 (t251 * t292 + (t264 * t249 + t269) * t248) * pkin(5) + t279, t279; 0, 0, 0, t254 * t238 + (-t258 * t239 - t291) * qJD(4) (-t239 * t272 + t248 * t275) * pkin(5) + t255, t255;];
JaD_transl  = t1;
