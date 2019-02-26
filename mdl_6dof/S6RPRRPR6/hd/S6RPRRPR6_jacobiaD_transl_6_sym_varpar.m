% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:48
% EndTime: 2019-02-26 21:03:48
% DurationCPUTime: 0.28s
% Computational Cost: add. (537->62), mult. (486->87), div. (0->0), fcn. (386->11), ass. (0->53)
t256 = qJ(4) + pkin(11);
t241 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t256);
t238 = t241 * qJD(4);
t242 = pkin(5) * cos(t256) + cos(qJ(4)) * pkin(4);
t240 = pkin(3) + t242;
t254 = pkin(10) + qJ(3);
t247 = sin(t254);
t249 = cos(t254);
t278 = t247 * qJD(5);
t291 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
t298 = t291 * t249;
t302 = (-t240 * t247 + t298) * qJD(3) + (t241 + pkin(7) + qJ(2)) * qJD(1) - t249 * t238 + t278;
t251 = qJ(6) + t256;
t244 = sin(t251);
t293 = r_i_i_C(2) * t244;
t245 = cos(t251);
t294 = r_i_i_C(1) * t245;
t267 = t240 - t293 + t294;
t264 = -t267 * t247 + t298;
t261 = cos(qJ(1));
t255 = qJD(4) + qJD(6);
t273 = t249 * t255 - qJD(1);
t300 = t261 * t273;
t292 = r_i_i_C(2) * t245;
t271 = r_i_i_C(1) * t244 + t292;
t297 = t271 * t255 + t238;
t285 = qJD(1) * t249;
t272 = -t255 + t285;
t259 = sin(qJ(1));
t280 = qJD(3) * t259;
t296 = -t247 * t280 + t272 * t261;
t289 = t247 * t255;
t279 = qJD(3) * t261;
t265 = t247 * t279 + t272 * t259;
t233 = t265 * t244 - t245 * t300;
t234 = t244 * t300 + t265 * t245;
t288 = t233 * r_i_i_C(1) + t234 * r_i_i_C(2);
t269 = t273 * t259;
t235 = t296 * t244 + t245 * t269;
t236 = t244 * t269 - t296 * t245;
t287 = -t235 * r_i_i_C(1) + t236 * r_i_i_C(2);
t284 = qJD(1) * t259;
t283 = qJD(1) * t261;
t282 = qJD(3) * t247;
t281 = qJD(3) * t249;
t276 = t291 * t247;
t270 = t241 * t285 - t238;
t239 = t242 * qJD(4);
t266 = qJD(1) * t242 - t239 * t249 + t241 * t282;
t263 = qJD(2) + t239 + (-t240 * t249 - cos(pkin(10)) * pkin(2) - pkin(1) - t276) * qJD(1);
t262 = qJD(5) * t249 + t297 * t247 + (-t267 * t249 - t276) * qJD(3);
t237 = t289 * t293;
t1 = [t236 * r_i_i_C(1) + t235 * r_i_i_C(2) - t302 * t259 + t263 * t261, t283, t262 * t261 - t264 * t284, t270 * t259 + t266 * t261 + t288, -t247 * t284 + t249 * t279, t288; -t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t263 * t259 + t302 * t261, t284, t262 * t259 + t264 * t283, t266 * t259 - t270 * t261 + t287, t247 * t283 + t249 * t280, t287; 0, 0, t264 * qJD(3) - t297 * t249 + t278, t237 + (-t255 * t294 - t239) * t247 + (-t241 - t271) * t281, t282, -t281 * t292 + t237 + (-t244 * t281 - t245 * t289) * r_i_i_C(1);];
JaD_transl  = t1;
