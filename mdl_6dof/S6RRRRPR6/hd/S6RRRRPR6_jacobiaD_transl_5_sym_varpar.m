% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:41
% EndTime: 2019-02-26 22:33:41
% DurationCPUTime: 0.27s
% Computational Cost: add. (462->65), mult. (522->93), div. (0->0), fcn. (410->10), ass. (0->58)
t246 = qJ(3) + qJ(4);
t241 = sin(t246);
t247 = sin(qJ(3));
t282 = pkin(3) * qJD(3);
t245 = qJD(3) + qJD(4);
t287 = pkin(4) * t245;
t231 = -t241 * t287 - t247 * t282;
t242 = cos(t246);
t239 = pkin(4) * t242;
t250 = cos(qJ(3));
t236 = t250 * pkin(3) + t239;
t234 = pkin(2) + t236;
t288 = pkin(4) * t241;
t235 = t247 * pkin(3) + t288;
t248 = sin(qJ(2));
t251 = cos(qJ(2));
t271 = t248 * qJD(5);
t283 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8);
t292 = t283 * t251;
t295 = (-t234 * t248 + t292) * qJD(2) + (pkin(7) + t235) * qJD(1) + t251 * t231 + t271;
t240 = pkin(11) + t246;
t237 = sin(t240);
t285 = r_i_i_C(2) * t237;
t238 = cos(t240);
t286 = r_i_i_C(1) * t238;
t259 = t234 - t285 + t286;
t254 = -t259 * t248 + t292;
t249 = sin(qJ(1));
t266 = t245 * t251 - qJD(1);
t261 = t266 * t249;
t252 = cos(qJ(1));
t262 = t266 * t252;
t264 = r_i_i_C(1) * t237 + r_i_i_C(2) * t238;
t291 = t264 * t245 - t231;
t276 = qJD(1) * t251;
t265 = -t245 + t276;
t274 = qJD(2) * t248;
t290 = -t249 * t274 + t265 * t252;
t280 = t245 * t248;
t272 = qJD(2) * t252;
t256 = t248 * t272 + t265 * t249;
t227 = t256 * t237 - t238 * t262;
t228 = t237 * t262 + t256 * t238;
t279 = t227 * r_i_i_C(1) + t228 * r_i_i_C(2);
t229 = t290 * t237 + t238 * t261;
t230 = t237 * t261 - t238 * t290;
t278 = -t229 * r_i_i_C(1) + t230 * r_i_i_C(2);
t277 = qJD(1) * t249;
t275 = qJD(1) * t252;
t273 = qJD(2) * t251;
t269 = t283 * t248;
t263 = t235 * t276 + t231;
t232 = t242 * t287 + t250 * t282;
t258 = qJD(1) * t236 - t232 * t251 + t235 * t274;
t255 = t232 + (-t234 * t251 - pkin(1) - t269) * qJD(1);
t253 = qJD(5) * t251 + t291 * t248 + (-t259 * t251 - t269) * qJD(2);
t233 = t280 * t285;
t1 = [t230 * r_i_i_C(1) + t229 * r_i_i_C(2) - t295 * t249 + t255 * t252, t253 * t252 - t254 * t277, t263 * t249 + t258 * t252 + t279 (t256 * t241 - t242 * t262) * pkin(4) + t279, -t248 * t277 + t251 * t272, 0; -t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t255 * t249 + t295 * t252, t253 * t249 + t254 * t275, t258 * t249 - t263 * t252 + t278 (-t241 * t290 - t242 * t261) * pkin(4) + t278, t248 * t275 + t249 * t273, 0; 0, t254 * qJD(2) - t291 * t251 + t271, t233 + (-t245 * t286 - t232) * t248 + (-t235 - t264) * t273, t233 + (-t286 - t239) * t280 + (-t264 - t288) * t273, t274, 0;];
JaD_transl  = t1;
