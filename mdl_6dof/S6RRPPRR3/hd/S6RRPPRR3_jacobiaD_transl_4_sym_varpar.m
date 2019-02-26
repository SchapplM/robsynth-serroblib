% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:34
% EndTime: 2019-02-26 21:29:34
% DurationCPUTime: 0.23s
% Computational Cost: add. (227->54), mult. (677->94), div. (0->0), fcn. (692->10), ass. (0->42)
t269 = sin(pkin(11));
t272 = cos(pkin(11));
t274 = sin(qJ(2));
t276 = cos(qJ(2));
t282 = t276 * t269 + t274 * t272;
t298 = qJD(2) * t282;
t268 = sin(pkin(12));
t270 = sin(pkin(6));
t271 = cos(pkin(12));
t273 = cos(pkin(6));
t297 = -t273 * t274 * pkin(2) + (r_i_i_C(1) * t268 + r_i_i_C(2) * t271 + pkin(8) + qJ(3)) * t270;
t296 = r_i_i_C(3) + qJ(4);
t275 = sin(qJ(1));
t295 = t274 * t275;
t277 = cos(qJ(1));
t294 = t274 * t277;
t293 = t275 * t276;
t292 = t276 * t277;
t291 = qJD(1) * t275;
t290 = qJD(2) * t274;
t289 = qJD(2) * t276;
t288 = pkin(2) * t290;
t287 = t273 * t289;
t261 = t274 * t269 - t276 * t272;
t256 = t261 * t273;
t285 = -t277 * t256 - t275 * t282;
t257 = t282 * t273;
t284 = -t277 * t257 + t275 * t261;
t283 = t275 * t257 + t277 * t261;
t281 = t271 * r_i_i_C(1) - t268 * r_i_i_C(2) + pkin(3);
t280 = t261 * qJD(2);
t253 = t273 * t269 * t290 - t272 * t287;
t259 = -t269 * t289 - t272 * t290;
t279 = t284 * qJD(1) + t275 * t253 + t277 * t259;
t278 = t283 * qJD(1) + t277 * t253 - t275 * t259;
t267 = t276 * pkin(2) + pkin(1);
t260 = pkin(2) * t287 - t270 * qJD(3);
t254 = t273 * t298;
t252 = t270 * t298;
t248 = t256 * t291 + (-qJD(1) * t282 - t254) * t277 + t275 * t280;
t245 = t285 * qJD(1) - t275 * t254 - t277 * t280;
t1 = [t285 * qJD(4) + t275 * t288 - t277 * t260 + t296 * t248 + t281 * t278 + (-t277 * t267 - t297 * t275) * qJD(1), -t283 * qJD(4) + t296 * t279 - t281 * t245 + ((t273 * t295 - t292) * qJD(2) + (-t273 * t292 + t295) * qJD(1)) * pkin(2), qJD(1) * t277 * t270, t245, 0, 0; -(t275 * t256 - t277 * t282) * qJD(4) - t277 * t288 - t275 * t260 + t296 * t245 + t281 * t279 + (-t275 * t267 + t297 * t277) * qJD(1), -t284 * qJD(4) - t296 * t278 + t281 * t248 + ((-t273 * t294 - t293) * qJD(2) + (-t273 * t293 - t294) * qJD(1)) * pkin(2), t270 * t291, -t248, 0, 0; 0, -t281 * t252 + (t282 * qJD(4) - t296 * t280 - t288) * t270, 0, t252, 0, 0;];
JaD_transl  = t1;
