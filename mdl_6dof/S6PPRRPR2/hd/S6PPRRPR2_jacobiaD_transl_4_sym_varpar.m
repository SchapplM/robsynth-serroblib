% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRPR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:57
% EndTime: 2019-02-26 19:40:58
% DurationCPUTime: 0.22s
% Computational Cost: add. (153->45), mult. (509->95), div. (0->0), fcn. (560->12), ass. (0->45)
t281 = sin(pkin(12));
t284 = sin(pkin(6));
t290 = sin(qJ(3));
t292 = cos(qJ(3));
t285 = cos(pkin(12));
t287 = cos(pkin(7));
t306 = t285 * t287;
t283 = sin(pkin(7));
t288 = cos(pkin(6));
t309 = t283 * t288;
t268 = (t281 * t292 + t290 * t306) * t284 + t290 * t309;
t286 = cos(pkin(11));
t282 = sin(pkin(11));
t311 = t282 * t288;
t277 = -t281 * t311 + t286 * t285;
t276 = -t286 * t281 - t285 * t311;
t310 = t283 * t284;
t296 = t276 * t287 + t282 * t310;
t264 = t277 * t292 + t296 * t290;
t305 = t286 * t288;
t275 = t281 * t305 + t282 * t285;
t274 = -t282 * t281 + t285 * t305;
t302 = t286 * t310;
t297 = -t274 * t287 + t302;
t316 = -t275 * t292 + t297 * t290;
t315 = -pkin(9) - r_i_i_C(3);
t314 = t275 * t290;
t308 = t284 * t285;
t307 = t284 * t287;
t304 = qJD(3) * t290;
t303 = qJD(3) * t292;
t300 = t283 * t303;
t299 = t287 * t303;
t289 = sin(qJ(4));
t291 = cos(qJ(4));
t298 = t289 * r_i_i_C(1) + t291 * r_i_i_C(2);
t294 = qJD(4) * t298;
t293 = qJD(3) * (t291 * r_i_i_C(1) - t289 * r_i_i_C(2) + pkin(3));
t273 = -t283 * t308 + t288 * t287;
t270 = -t276 * t283 + t282 * t307;
t269 = -t274 * t283 - t286 * t307;
t265 = t284 * t281 * t304 - t288 * t300 - t299 * t308;
t259 = -t282 * t284 * t300 - t276 * t299 + t277 * t304;
t257 = -t274 * t299 + (t292 * t302 + t314) * qJD(3);
t1 = [0, 0, t315 * t259 - (-t277 * t290 + t296 * t292) * t294 - t264 * t293, t298 * t259 + ((-t264 * t291 - t270 * t289) * r_i_i_C(1) + (t264 * t289 - t270 * t291) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, t315 * t257 - (-t297 * t292 - t314) * t294 + t316 * t293, t298 * t257 + ((-t269 * t289 + t291 * t316) * r_i_i_C(1) + (-t269 * t291 - t289 * t316) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, t315 * t265 - (t292 * t309 + (-t281 * t290 + t292 * t306) * t284) * t294 - t268 * t293, t298 * t265 + ((-t268 * t291 - t273 * t289) * r_i_i_C(1) + (t268 * t289 - t273 * t291) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
