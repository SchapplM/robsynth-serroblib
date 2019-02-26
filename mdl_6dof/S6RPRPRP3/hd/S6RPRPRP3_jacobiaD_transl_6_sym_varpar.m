% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP3
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

function JaD_transl = S6RPRPRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:52
% EndTime: 2019-02-26 20:44:53
% DurationCPUTime: 0.36s
% Computational Cost: add. (499->64), mult. (547->96), div. (0->0), fcn. (460->10), ass. (0->46)
t266 = sin(qJ(3));
t257 = cos(pkin(10)) * pkin(4) + pkin(3);
t262 = pkin(10) + qJ(5);
t258 = sin(t262);
t260 = cos(t262);
t298 = r_i_i_C(3) + qJ(6);
t300 = r_i_i_C(1) + pkin(5);
t301 = t298 * t258 + t300 * t260 + t257;
t267 = cos(qJ(3));
t299 = r_i_i_C(2) + pkin(8) + qJ(4);
t303 = t299 * t267;
t307 = t301 * t266 - t303;
t284 = t266 * qJD(4);
t286 = qJD(6) * t258;
t306 = (-t257 * t266 + t303) * qJD(3) + t267 * t286 + t284;
t304 = -t286 + (t300 * t258 - t298 * t260) * qJD(5);
t263 = qJ(1) + pkin(9);
t259 = sin(t263);
t296 = t259 * t267;
t261 = cos(t263);
t295 = t261 * t260;
t294 = qJD(1) * t259;
t293 = qJD(1) * t261;
t292 = qJD(1) * t266;
t291 = qJD(3) * t266;
t290 = qJD(3) * t267;
t289 = qJD(5) * t259;
t288 = qJD(5) * t261;
t287 = qJD(5) * t266;
t285 = qJD(6) * t260;
t283 = pkin(4) * sin(pkin(10)) + pkin(7);
t282 = t259 * t291;
t281 = t261 * t291;
t280 = t258 * t289;
t279 = t258 * t288;
t278 = t260 * t288;
t277 = t299 * t266;
t274 = t259 * t258 + t267 * t295;
t272 = -t257 * t267 - pkin(2) - t277;
t271 = t258 * t293 + t260 * t289;
t268 = qJD(4) * t267 + t304 * t266 + (-t267 * t301 - t277) * qJD(3);
t246 = t274 * qJD(1) - t260 * t282 - t267 * t280 - t278;
t245 = -t258 * t282 - t260 * t294 + t271 * t267 - t279;
t244 = t267 * t279 + (t267 * t294 + t281) * t260 - t271;
t243 = t258 * t281 - t267 * t278 - t280 + (t258 * t296 + t295) * qJD(1);
t1 = [-t261 * t285 - t300 * t246 - t298 * t245 - t306 * t259 + (-cos(qJ(1)) * pkin(1) - t283 * t259 + t272 * t261) * qJD(1), 0, t268 * t261 + t307 * t294, -t259 * t292 + t261 * t290, t274 * qJD(6) + t300 * t243 - t298 * t244, -t243; -t259 * t285 - t300 * t244 - t298 * t243 + t306 * t261 + (-sin(qJ(1)) * pkin(1) + t283 * t261 + t272 * t259) * qJD(1), 0, t268 * t259 - t293 * t307, t259 * t290 + t261 * t292 -(t261 * t258 - t260 * t296) * qJD(6) + t298 * t246 - t300 * t245, t245; 0, 0, -qJD(3) * t307 - t304 * t267 + t284, t291 (-t298 * t287 - t300 * t290) * t258 + (t298 * t290 + (-t300 * qJD(5) + qJD(6)) * t266) * t260, t258 * t290 + t260 * t287;];
JaD_transl  = t1;
