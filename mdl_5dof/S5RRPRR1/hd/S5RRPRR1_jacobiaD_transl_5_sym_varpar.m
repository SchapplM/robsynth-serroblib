% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:53
% EndTime: 2019-07-18 17:22:53
% DurationCPUTime: 0.27s
% Computational Cost: add. (243->54), mult. (358->88), div. (0->0), fcn. (275->8), ass. (0->51)
t255 = qJ(2) + qJ(4);
t252 = sin(t255);
t257 = sin(qJ(5));
t260 = cos(qJ(5));
t289 = qJD(5) * t260;
t253 = cos(t255);
t254 = qJD(2) + qJD(4);
t297 = t253 * t254;
t304 = t252 * t289 + t257 * t297;
t300 = pkin(4) + r_i_i_C(3);
t303 = t300 * t297;
t288 = r_i_i_C(2) * t252 * t257;
t302 = t300 * t253 + t288;
t290 = qJD(5) * t257;
t279 = t252 * t290;
t301 = r_i_i_C(1) * t279 + t304 * r_i_i_C(2);
t299 = r_i_i_C(1) * t260;
t259 = sin(qJ(1));
t298 = t252 * t259;
t296 = t254 * t260;
t258 = sin(qJ(2));
t263 = pkin(2) + pkin(1);
t295 = t258 * t263;
t262 = cos(qJ(1));
t294 = t260 * t262;
t293 = qJD(1) * t259;
t292 = qJD(1) * t262;
t291 = qJD(2) * t263;
t287 = qJD(1) * t299;
t286 = t300 * t252;
t285 = t300 * t254;
t284 = t252 * t296;
t282 = t253 * t296;
t280 = t258 * t291;
t274 = qJD(5) * t253 - qJD(1);
t273 = qJD(1) * t253 - qJD(5);
t272 = t301 * t262 + t287 * t298;
t271 = t274 * t257;
t270 = t301 * t259 + t302 * t292;
t261 = cos(qJ(2));
t269 = -t261 * t263 - t286;
t267 = (-t253 * t299 - t286) * t254;
t266 = t252 * t254 * t262 + t273 * t259;
t265 = -t261 * t291 + t267;
t264 = -t253 * r_i_i_C(2) * t289 + t254 * t288 + (-t253 * t290 - t284) * r_i_i_C(1) + t303;
t256 = -pkin(3) - qJ(3);
t238 = -t273 * t294 + (t271 + t284) * t259;
t237 = t274 * t260 * t259 + (-t254 * t298 + t273 * t262) * t257;
t236 = t266 * t260 + t262 * t271;
t235 = t266 * t257 - t274 * t294;
t1 = [t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t262 * qJD(3) + (-t253 * t285 + t280) * t259 + (t256 * t259 + t269 * t262) * qJD(1), t265 * t262 + (-t302 + t295) * t293 + t272, t292, t262 * t267 - t293 * t302 + t272, t235 * r_i_i_C(1) + t236 * r_i_i_C(2); -t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t259 * qJD(3) + (-t280 + t303) * t262 + (-t256 * t262 + t269 * t259) * qJD(1), (-t252 * t299 - t295) * t292 + t265 * t259 + t270, t293, -t259 * r_i_i_C(1) * t282 + (-t259 * t285 - t262 * t287) * t252 + t270, -t237 * r_i_i_C(1) + t238 * r_i_i_C(2); 0, t264 - t280, 0, t264, (t279 - t282) * r_i_i_C(2) - t304 * r_i_i_C(1);];
JaD_transl  = t1;
