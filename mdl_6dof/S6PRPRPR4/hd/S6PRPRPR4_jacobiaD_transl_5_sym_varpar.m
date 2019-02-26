% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:17
% EndTime: 2019-02-26 19:48:17
% DurationCPUTime: 0.18s
% Computational Cost: add. (258->45), mult. (529->81), div. (0->0), fcn. (515->11), ass. (0->40)
t274 = sin(pkin(10));
t277 = cos(pkin(10));
t281 = cos(qJ(2));
t278 = cos(pkin(6));
t280 = sin(qJ(2));
t296 = t278 * t280;
t263 = t274 * t281 + t277 * t296;
t272 = pkin(11) + qJ(4);
t270 = sin(t272);
t271 = cos(t272);
t275 = sin(pkin(6));
t298 = t275 * t277;
t302 = -t263 * t271 + t270 * t298;
t301 = r_i_i_C(3) + qJ(5);
t299 = t274 * t275;
t297 = t275 * t280;
t295 = t278 * t281;
t294 = qJD(2) * t280;
t293 = qJD(2) * t281;
t291 = t275 * t293;
t290 = t274 * t294;
t289 = t277 * t293;
t273 = sin(pkin(12));
t276 = cos(pkin(12));
t288 = -t276 * r_i_i_C(1) + t273 * r_i_i_C(2) - pkin(4);
t287 = t273 * r_i_i_C(1) + t276 * r_i_i_C(2) + pkin(8) + qJ(3);
t265 = -t274 * t296 + t277 * t281;
t286 = t265 * t271 + t270 * t299;
t285 = t278 * t270 + t271 * t297;
t284 = t274 * t295 + t277 * t280;
t283 = -t301 * t270 + t288 * t271 - cos(pkin(11)) * pkin(3) - pkin(2);
t282 = t270 * qJD(5) + (t288 * t270 + t301 * t271) * qJD(4);
t261 = -t278 * t290 + t289;
t260 = t284 * qJD(2);
t259 = t263 * qJD(2);
t258 = -t278 * t289 + t290;
t256 = t285 * qJD(4) + t270 * t291;
t254 = t286 * qJD(4) - t260 * t270;
t252 = -t302 * qJD(4) - t258 * t270;
t1 = [0, t265 * qJD(3) - t287 * t260 + t283 * t261 - t282 * t284, t261, t286 * qJD(5) + t301 * (-t260 * t271 + (-t265 * t270 + t271 * t299) * qJD(4)) + t288 * t254, t254, 0; 0, t263 * qJD(3) - t287 * t258 + t282 * (-t274 * t280 + t277 * t295) + t283 * t259, t259, -t302 * qJD(5) + t301 * (-t258 * t271 + (-t263 * t270 - t271 * t298) * qJD(4)) + t288 * t252, t252, 0; 0 (t280 * qJD(3) + t282 * t281 + (t283 * t280 + t287 * t281) * qJD(2)) * t275, t275 * t294, t285 * qJD(5) + t301 * (t271 * t291 + (-t270 * t297 + t271 * t278) * qJD(4)) + t288 * t256, t256, 0;];
JaD_transl  = t1;
