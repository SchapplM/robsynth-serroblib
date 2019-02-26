% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:52
% EndTime: 2019-02-26 21:38:52
% DurationCPUTime: 0.33s
% Computational Cost: add. (544->63), mult. (506->88), div. (0->0), fcn. (399->12), ass. (0->53)
t259 = qJ(4) + pkin(11);
t245 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t259);
t242 = t245 * qJD(4);
t246 = pkin(5) * cos(t259) + cos(qJ(4)) * pkin(4);
t244 = pkin(3) + t246;
t260 = qJ(2) + pkin(10);
t252 = sin(t260);
t254 = cos(t260);
t297 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
t274 = t297 * t254 - sin(qJ(2)) * pkin(2);
t284 = t252 * qJD(5);
t310 = (-t244 * t252 + t274) * qJD(2) + (t245 + qJ(3) + pkin(7)) * qJD(1) - t254 * t242 + t284;
t255 = qJ(6) + t259;
t248 = sin(t255);
t300 = r_i_i_C(2) * t248;
t249 = cos(t255);
t301 = r_i_i_C(1) * t249;
t275 = t244 - t300 + t301;
t269 = -t275 * t252 + t274;
t267 = cos(qJ(1));
t258 = qJD(4) + qJD(6);
t281 = t254 * t258 - qJD(1);
t308 = t267 * t281;
t306 = -t297 * t252 - cos(qJ(2)) * pkin(2);
t299 = r_i_i_C(2) * t249;
t279 = r_i_i_C(1) * t248 + t299;
t305 = t279 * t258 + t242;
t291 = qJD(1) * t254;
t280 = -t258 + t291;
t264 = sin(qJ(1));
t286 = qJD(2) * t264;
t304 = -t252 * t286 + t280 * t267;
t295 = t252 * t258;
t285 = qJD(2) * t267;
t271 = t252 * t285 + t280 * t264;
t237 = t271 * t248 - t249 * t308;
t238 = t248 * t308 + t271 * t249;
t294 = t237 * r_i_i_C(1) + t238 * r_i_i_C(2);
t277 = t281 * t264;
t239 = t304 * t248 + t249 * t277;
t240 = t248 * t277 - t304 * t249;
t293 = -t239 * r_i_i_C(1) + t240 * r_i_i_C(2);
t290 = qJD(1) * t264;
t289 = qJD(1) * t267;
t288 = qJD(2) * t252;
t287 = qJD(2) * t254;
t278 = t245 * t291 - t242;
t243 = t246 * qJD(4);
t272 = qJD(1) * t246 - t243 * t254 + t245 * t288;
t270 = qJD(3) + t243 + (-t244 * t254 - pkin(1) + t306) * qJD(1);
t268 = qJD(5) * t254 + t305 * t252 + (-t275 * t254 + t306) * qJD(2);
t241 = t295 * t300;
t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) - t264 * t310 + t270 * t267, t268 * t267 - t269 * t290, t289, t278 * t264 + t272 * t267 + t294, -t252 * t290 + t254 * t285, t294; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t270 * t264 + t267 * t310, t268 * t264 + t269 * t289, t290, t272 * t264 - t278 * t267 + t293, t252 * t289 + t254 * t286, t293; 0, t269 * qJD(2) - t305 * t254 + t284, 0, t241 + (-t258 * t301 - t243) * t252 + (-t245 - t279) * t287, t288, -t287 * t299 + t241 + (-t248 * t287 - t249 * t295) * r_i_i_C(1);];
JaD_transl  = t1;
