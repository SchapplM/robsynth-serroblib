% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:51
% EndTime: 2019-02-26 21:24:52
% DurationCPUTime: 0.30s
% Computational Cost: add. (489->62), mult. (569->89), div. (0->0), fcn. (477->10), ass. (0->46)
t267 = cos(pkin(10)) * pkin(4) + pkin(3);
t274 = qJ(2) + pkin(9);
t270 = sin(t274);
t272 = cos(t274);
t314 = r_i_i_C(2) + pkin(8) + qJ(4);
t292 = t314 * t272 - sin(qJ(2)) * pkin(2);
t301 = t270 * qJD(4);
t273 = pkin(10) + qJ(5);
t269 = sin(t273);
t302 = qJD(6) * t269;
t324 = (-t267 * t270 + t292) * qJD(2) + (pkin(4) * sin(pkin(10)) + qJ(3) + pkin(7)) * qJD(1) + t272 * t302 + t301;
t271 = cos(t273);
t313 = r_i_i_C(3) + qJ(6);
t317 = pkin(5) + r_i_i_C(1);
t288 = -t313 * t269 - t317 * t271;
t285 = -t267 + t288;
t283 = t285 * t270 + t292;
t318 = t317 * t269 - t313 * t271;
t322 = t318 * qJD(5) - t302;
t320 = -t314 * t270 - cos(qJ(2)) * pkin(2);
t281 = cos(qJ(1));
t311 = t271 * t281;
t279 = sin(qJ(1));
t310 = t279 * t269;
t309 = qJD(1) * t279;
t308 = qJD(1) * t281;
t307 = qJD(2) * t272;
t306 = qJD(2) * t279;
t305 = qJD(2) * t281;
t304 = qJD(5) * t279;
t303 = qJD(5) * t281;
t300 = t271 * qJD(6);
t299 = t270 * t306;
t298 = t270 * t305;
t297 = t269 * t304;
t296 = t269 * t303;
t295 = t271 * t303;
t290 = t272 * t311 + t310;
t286 = t269 * t308 + t271 * t304;
t284 = -t300 + qJD(3) + (-t267 * t272 - pkin(1) + t320) * qJD(1);
t282 = qJD(4) * t272 + t322 * t270 + (t285 * t272 + t320) * qJD(2);
t256 = t290 * qJD(1) - t271 * t299 - t272 * t297 - t295;
t255 = -t269 * t299 - t271 * t309 + t286 * t272 - t296;
t254 = t272 * t296 + (t272 * t309 + t298) * t271 - t286;
t253 = t269 * t298 - t272 * t295 - t297 + (t272 * t310 + t311) * qJD(1);
t1 = [-t313 * t255 - t317 * t256 - t324 * t279 + t284 * t281, t282 * t281 - t283 * t309, t308, -t270 * t309 + t272 * t305, t290 * qJD(6) + t317 * t253 - t313 * t254, -t253; -t313 * t253 - t317 * t254 + t284 * t279 + t324 * t281, t282 * t279 + t283 * t308, t309, t270 * t308 + t272 * t306 -(-t279 * t272 * t271 + t269 * t281) * qJD(6) + t313 * t256 - t317 * t255, t255; 0, t283 * qJD(2) - t322 * t272 + t301, 0, qJD(2) * t270, -t318 * t307 + (t288 * qJD(5) + t300) * t270, qJD(5) * t270 * t271 + t269 * t307;];
JaD_transl  = t1;
