% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP2
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
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:16
% EndTime: 2019-02-26 20:44:17
% DurationCPUTime: 0.30s
% Computational Cost: add. (472->59), mult. (540->85), div. (0->0), fcn. (448->10), ass. (0->46)
t272 = qJ(3) + pkin(10);
t268 = sin(t272);
t270 = cos(t272);
t312 = pkin(8) + r_i_i_C(2);
t290 = t312 * t270 - sin(qJ(3)) * pkin(3);
t318 = -pkin(4) * t268 + t290;
t275 = sin(qJ(5));
t277 = cos(qJ(5));
t308 = r_i_i_C(3) + qJ(6);
t313 = pkin(5) + r_i_i_C(1);
t285 = -t308 * t275 - t313 * t277;
t281 = -pkin(4) + t285;
t280 = t281 * t268 + t290;
t314 = t313 * t275 - t308 * t277;
t317 = t314 * qJD(5) - qJD(6) * t275;
t315 = -t312 * t268 - cos(qJ(3)) * pkin(3);
t273 = qJ(1) + pkin(9);
t269 = sin(t273);
t307 = t269 * t275;
t306 = t269 * t277;
t271 = cos(t273);
t305 = t271 * t275;
t304 = t271 * t277;
t303 = qJD(1) * t269;
t302 = qJD(1) * t271;
t301 = qJD(3) * t269;
t300 = qJD(3) * t270;
t299 = qJD(3) * t271;
t298 = qJD(5) * t275;
t297 = qJD(5) * t277;
t294 = t268 * t299;
t293 = t268 * t301;
t292 = t269 * t298;
t291 = t271 * t297;
t288 = t270 * t304 + t307;
t287 = t270 * t307 + t304;
t286 = -pkin(4) * t270 - pkin(2) + t315;
t283 = t269 * t297 + t275 * t302;
t282 = t271 * t298 + t277 * t303;
t279 = t317 * t268 + (t281 * t270 + t315) * qJD(3);
t274 = -qJ(4) - pkin(7);
t256 = t288 * qJD(1) - t270 * t292 - t277 * t293 - t291;
t255 = t283 * t270 - t275 * t293 - t282;
t254 = t282 * t270 + t277 * t294 - t283;
t253 = t287 * qJD(1) - t270 * t291 + t275 * t294 - t292;
t1 = [-t287 * qJD(6) + t271 * qJD(4) - t313 * t256 - t308 * t255 - t318 * t301 + (-cos(qJ(1)) * pkin(1) + t269 * t274 + t286 * t271) * qJD(1), 0, t279 * t271 - t280 * t303, t302, t288 * qJD(6) + t313 * t253 - t308 * t254, -t253; -(-t270 * t305 + t306) * qJD(6) + t269 * qJD(4) - t313 * t254 - t308 * t253 + t318 * t299 + (-sin(qJ(1)) * pkin(1) - t271 * t274 + t286 * t269) * qJD(1), 0, t279 * t269 + t280 * t302, t303 -(-t270 * t306 + t305) * qJD(6) + t308 * t256 - t313 * t255, t255; 0, 0, t280 * qJD(3) - t317 * t270, 0, -t314 * t300 + (t285 * qJD(5) + t277 * qJD(6)) * t268, t268 * t297 + t275 * t300;];
JaD_transl  = t1;
