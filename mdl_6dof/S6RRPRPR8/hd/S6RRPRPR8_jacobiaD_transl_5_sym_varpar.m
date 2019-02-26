% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:48
% EndTime: 2019-02-26 21:41:49
% DurationCPUTime: 0.31s
% Computational Cost: add. (335->61), mult. (543->94), div. (0->0), fcn. (458->8), ass. (0->45)
t259 = sin(qJ(2));
t253 = cos(pkin(10)) * pkin(3) + pkin(2);
t256 = pkin(10) + qJ(4);
t254 = sin(t256);
t255 = cos(t256);
t293 = r_i_i_C(3) + qJ(5);
t295 = pkin(4) + r_i_i_C(1);
t296 = t293 * t254 + t295 * t255 + t253;
t261 = cos(qJ(2));
t294 = r_i_i_C(2) + pkin(8) + qJ(3);
t298 = t294 * t261;
t302 = t296 * t259 - t298;
t279 = t259 * qJD(3);
t281 = qJD(5) * t254;
t301 = (-t253 * t259 + t298) * qJD(2) + t261 * t281 + t279;
t299 = -t281 + (t295 * t254 - t293 * t255) * qJD(4);
t260 = sin(qJ(1));
t291 = t260 * t261;
t262 = cos(qJ(1));
t290 = t262 * t255;
t289 = qJD(1) * t260;
t288 = qJD(1) * t262;
t287 = qJD(2) * t259;
t286 = qJD(2) * t261;
t285 = qJD(2) * t262;
t284 = qJD(4) * t259;
t283 = qJD(4) * t260;
t282 = qJD(4) * t262;
t280 = t255 * qJD(5);
t278 = pkin(3) * sin(pkin(10)) + pkin(7);
t277 = t260 * t287;
t276 = t254 * t283;
t275 = t259 * t285;
t274 = t254 * t282;
t273 = t255 * t282;
t272 = t294 * t259;
t269 = t260 * t254 + t261 * t290;
t267 = -t253 * t261 - pkin(1) - t272;
t266 = t254 * t288 + t255 * t283;
t263 = qJD(3) * t261 + t299 * t259 + (-t261 * t296 - t272) * qJD(2);
t242 = qJD(1) * t269 - t255 * t277 - t261 * t276 - t273;
t241 = -t254 * t277 - t255 * t289 + t261 * t266 - t274;
t240 = t261 * t274 + (t261 * t289 + t275) * t255 - t266;
t239 = t254 * t275 - t261 * t273 - t276 + (t254 * t291 + t290) * qJD(1);
t1 = [-t262 * t280 - t295 * t242 - t293 * t241 - t301 * t260 + (-t260 * t278 + t262 * t267) * qJD(1), t263 * t262 + t302 * t289, -t259 * t289 + t261 * t285, t269 * qJD(5) + t295 * t239 - t293 * t240, -t239, 0; -t260 * t280 - t295 * t240 - t293 * t239 + t301 * t262 + (t260 * t267 + t262 * t278) * qJD(1), t260 * t263 - t288 * t302, t259 * t288 + t260 * t286 -(t262 * t254 - t255 * t291) * qJD(5) + t293 * t242 - t295 * t241, t241, 0; 0, -qJD(2) * t302 - t299 * t261 + t279, t287 (-t293 * t284 - t295 * t286) * t254 + (t293 * t286 + (-t295 * qJD(4) + qJD(5)) * t259) * t255, t254 * t286 + t255 * t284, 0;];
JaD_transl  = t1;
