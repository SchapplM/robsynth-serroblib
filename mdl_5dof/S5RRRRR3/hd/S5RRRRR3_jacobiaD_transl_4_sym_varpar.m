% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:39
% EndTime: 2019-07-18 17:19:39
% DurationCPUTime: 0.26s
% Computational Cost: add. (264->54), mult. (390->92), div. (0->0), fcn. (295->8), ass. (0->53)
t256 = qJ(2) + qJ(3);
t253 = sin(t256);
t260 = cos(qJ(4));
t307 = r_i_i_C(1) * t260 + pkin(2);
t310 = t253 * t307;
t257 = sin(qJ(4));
t290 = qJD(4) * t260;
t254 = cos(t256);
t255 = qJD(2) + qJD(3);
t298 = t254 * t255;
t309 = t253 * t290 + t257 * t298;
t304 = pkin(5) + r_i_i_C(3);
t285 = t304 * t254;
t258 = sin(qJ(2));
t299 = pkin(1) * qJD(2);
t287 = t258 * t299;
t302 = pkin(2) * t253;
t308 = -t287 + (t285 - t302) * t255;
t291 = qJD(4) * t257;
t278 = t253 * t291;
t305 = r_i_i_C(1) * t278 + t309 * r_i_i_C(2);
t303 = pkin(1) * t258;
t300 = r_i_i_C(2) * t257;
t259 = sin(qJ(1));
t297 = t255 * t259;
t296 = t255 * t260;
t262 = cos(qJ(1));
t295 = t255 * t262;
t294 = t260 * t262;
t293 = qJD(1) * t259;
t292 = qJD(1) * t262;
t289 = t253 * t300;
t288 = qJD(1) * t300;
t286 = t304 * t253;
t284 = t304 * t259;
t283 = t253 * t296;
t273 = qJD(4) * t254 - qJD(1);
t272 = qJD(1) * t254 - qJD(4);
t271 = t307 * t254;
t270 = t307 * t262;
t269 = t305 * t262 + t293 * t310;
t268 = t273 * t257;
t267 = t262 * t253 * t288 + t305 * t259 + t292 * t285;
t266 = t253 * t295 + t272 * t259;
t261 = cos(qJ(2));
t265 = qJD(1) * (-pkin(1) * t261 - pkin(2) * t254 - t286);
t264 = -t261 * t299 + (-t271 - t286) * t255;
t263 = -t254 * r_i_i_C(2) * t290 + (-t254 * t291 - t283) * r_i_i_C(1) + t304 * t298 + (-t302 + t289) * t255;
t238 = -t272 * t294 + (t268 + t283) * t259;
t237 = t273 * t260 * t259 + (-t253 * t297 + t272 * t262) * t257;
t236 = t266 * t260 + t262 * t268;
t235 = t266 * t257 - t273 * t294;
t1 = [t238 * r_i_i_C(1) + t237 * r_i_i_C(2) - t308 * t259 + t262 * t265, (-t285 - t289 + t303) * t293 + t264 * t262 + t269, (-t259 * t288 - t304 * t295) * t253 + (-qJD(1) * t284 - t255 * t270) * t254 + t269, t235 * r_i_i_C(1) + t236 * r_i_i_C(2), 0; -t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t259 * t265 + t308 * t262, (-t303 - t310) * t292 + t264 * t259 + t267, -t271 * t297 + (-qJD(1) * t270 - t255 * t284) * t253 + t267, -t237 * r_i_i_C(1) + t238 * r_i_i_C(2), 0; 0, t263 - t287, t263, (-t254 * t296 + t278) * r_i_i_C(2) - t309 * r_i_i_C(1), 0;];
JaD_transl  = t1;
