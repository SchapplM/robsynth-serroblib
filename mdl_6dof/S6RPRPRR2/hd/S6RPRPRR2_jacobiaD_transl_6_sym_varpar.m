% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:39
% EndTime: 2019-02-26 20:49:40
% DurationCPUTime: 0.29s
% Computational Cost: add. (504->60), mult. (450->90), div. (0->0), fcn. (353->12), ass. (0->52)
t267 = cos(qJ(5));
t252 = pkin(5) * t267 + pkin(4);
t261 = qJ(3) + pkin(11);
t254 = sin(t261);
t256 = cos(t261);
t300 = r_i_i_C(3) + pkin(9) + pkin(8);
t276 = t300 * t256 - sin(qJ(3)) * pkin(3);
t265 = sin(qJ(5));
t299 = pkin(5) * qJD(5);
t289 = t265 * t299;
t312 = (-t252 * t254 + t276) * qJD(3) - t256 * t289;
t263 = qJ(5) + qJ(6);
t259 = cos(t263);
t258 = sin(t263);
t301 = r_i_i_C(2) * t258;
t277 = r_i_i_C(1) * t259 + t252 - t301;
t271 = -t277 * t254 + t276;
t260 = qJD(5) + qJD(6);
t293 = qJD(1) * t256;
t283 = -t260 + t293;
t310 = t259 * t283;
t309 = t265 * (-qJD(5) + t293);
t307 = -t300 * t254 - cos(qJ(3)) * pkin(3);
t284 = t256 * t260 - qJD(1);
t291 = qJD(3) * t258;
t306 = -t254 * t291 + t284 * t259;
t279 = r_i_i_C(1) * t258 + r_i_i_C(2) * t259;
t305 = t279 * t260 + t289;
t302 = pkin(5) * t265;
t297 = t259 * t260;
t262 = qJ(1) + pkin(10);
t255 = sin(t262);
t257 = cos(t262);
t278 = t283 * t258;
t247 = t255 * t278 - t306 * t257;
t290 = qJD(3) * t259;
t273 = t254 * t290 + t284 * t258;
t248 = t255 * t310 + t273 * t257;
t296 = t247 * r_i_i_C(1) + t248 * r_i_i_C(2);
t249 = t306 * t255 + t257 * t278;
t250 = t273 * t255 - t257 * t310;
t295 = -t249 * r_i_i_C(1) + t250 * r_i_i_C(2);
t294 = qJD(1) * t255;
t292 = qJD(1) * t257;
t288 = t267 * t299;
t285 = qJ(4) + pkin(7) + t302;
t280 = qJD(4) + t288;
t274 = -t252 * t256 - pkin(2) + t307;
t272 = qJD(3) * t254 * t265 + (-qJD(5) * t256 + qJD(1)) * t267;
t270 = t305 * t254 + (-t277 * t256 + t307) * qJD(3);
t251 = t254 * t260 * t301;
t1 = [t250 * r_i_i_C(1) + t249 * r_i_i_C(2) + t280 * t257 - t312 * t255 + (-cos(qJ(1)) * pkin(1) - t285 * t255 + t274 * t257) * qJD(1), 0, t270 * t257 - t271 * t294, t292 (t255 * t309 + t272 * t257) * pkin(5) + t296, t296; -t248 * r_i_i_C(1) + t247 * r_i_i_C(2) + t280 * t255 + t312 * t257 + (-sin(qJ(1)) * pkin(1) + t285 * t257 + t274 * t255) * qJD(1), 0, t270 * t255 + t271 * t292, t294 (t272 * t255 - t257 * t309) * pkin(5) + t295, t295; 0, 0, t271 * qJD(3) - t305 * t256, 0, t251 + (-r_i_i_C(1) * t297 - t288) * t254 + (-t279 - t302) * t256 * qJD(3), -t256 * r_i_i_C(2) * t290 + t251 + (-t254 * t297 - t256 * t291) * r_i_i_C(1);];
JaD_transl  = t1;
