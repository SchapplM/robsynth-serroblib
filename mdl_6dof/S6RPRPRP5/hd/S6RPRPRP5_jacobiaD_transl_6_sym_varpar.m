% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:54
% EndTime: 2019-02-26 20:45:54
% DurationCPUTime: 0.28s
% Computational Cost: add. (482->63), mult. (549->92), div. (0->0), fcn. (464->9), ass. (0->48)
t264 = cos(pkin(10)) * pkin(4) + pkin(3);
t271 = pkin(9) + qJ(3);
t267 = sin(t271);
t269 = cos(t271);
t297 = t267 * qJD(4);
t270 = pkin(10) + qJ(5);
t266 = sin(t270);
t298 = qJD(6) * t266;
t310 = r_i_i_C(2) + pkin(8) + qJ(4);
t313 = t310 * t269;
t317 = (-t264 * t267 + t313) * qJD(3) + t269 * t298 + t297;
t268 = cos(t270);
t309 = r_i_i_C(3) + qJ(6);
t311 = pkin(5) + r_i_i_C(1);
t312 = t266 * t311 - t268 * t309;
t316 = qJD(5) * t312 - t298;
t283 = -t266 * t309 - t268 * t311;
t279 = -t264 + t283;
t278 = t279 * t267 + t313;
t275 = sin(qJ(1));
t307 = t275 * t266;
t276 = cos(qJ(1));
t306 = t276 * t268;
t305 = qJD(1) * t275;
t304 = qJD(1) * t276;
t303 = qJD(3) * t269;
t302 = qJD(3) * t275;
t301 = qJD(3) * t276;
t300 = qJD(5) * t275;
t299 = qJD(5) * t276;
t296 = t268 * qJD(6);
t295 = t267 * t302;
t294 = t266 * t300;
t293 = t267 * t301;
t292 = t266 * t299;
t291 = t268 * t299;
t290 = t310 * t267;
t287 = pkin(4) * sin(pkin(10)) + pkin(7) + qJ(2);
t286 = qJD(2) - t296;
t285 = t269 * t306 + t307;
t281 = -t264 * t269 - cos(pkin(9)) * pkin(2) - pkin(1) - t290;
t280 = t266 * t304 + t268 * t300;
t277 = qJD(4) * t269 + t316 * t267 + (t269 * t279 - t290) * qJD(3);
t253 = qJD(1) * t285 - t268 * t295 - t269 * t294 - t291;
t252 = -t266 * t295 - t268 * t305 + t269 * t280 - t292;
t251 = t269 * t292 + (t269 * t305 + t293) * t268 - t280;
t250 = t266 * t293 - t269 * t291 - t294 + (t269 * t307 + t306) * qJD(1);
t1 = [t286 * t276 - t311 * t253 - t309 * t252 - t317 * t275 + (-t275 * t287 + t276 * t281) * qJD(1), t304, t276 * t277 - t278 * t305, -t267 * t305 + t269 * t301, qJD(6) * t285 + t250 * t311 - t251 * t309, -t250; t286 * t275 - t311 * t251 - t309 * t250 + t317 * t276 + (t275 * t281 + t276 * t287) * qJD(1), t305, t275 * t277 + t278 * t304, t267 * t304 + t269 * t302 -(-t268 * t269 * t275 + t266 * t276) * qJD(6) + t309 * t253 - t311 * t252, t252; 0, 0, qJD(3) * t278 - t269 * t316 + t297, qJD(3) * t267, -t312 * t303 + (qJD(5) * t283 + t296) * t267, qJD(5) * t267 * t268 + t266 * t303;];
JaD_transl  = t1;
