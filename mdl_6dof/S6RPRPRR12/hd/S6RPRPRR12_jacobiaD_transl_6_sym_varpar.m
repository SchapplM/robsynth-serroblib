% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:03
% EndTime: 2019-02-26 20:55:03
% DurationCPUTime: 0.26s
% Computational Cost: add. (283->61), mult. (498->92), div. (0->0), fcn. (390->8), ass. (0->52)
t232 = qJD(5) + qJD(6);
t237 = cos(qJ(5));
t268 = pkin(5) * qJD(5);
t249 = t237 * t268 + qJD(4);
t233 = qJ(5) + qJ(6);
t230 = sin(t233);
t270 = r_i_i_C(2) * t230;
t231 = cos(t233);
t271 = r_i_i_C(1) * t231;
t242 = (-t270 + t271) * t232 + t249;
t259 = pkin(3) + r_i_i_C(3) + pkin(9) + pkin(8);
t280 = -t259 * qJD(3) + t242;
t235 = sin(qJ(3));
t238 = cos(qJ(3));
t234 = sin(qJ(5));
t255 = pkin(5) * t234 + qJ(4);
t258 = t234 * t268;
t279 = -t258 + (t259 * t235 - t255 * t238 + qJ(2)) * qJD(1);
t239 = cos(qJ(1));
t254 = t232 * t238 + qJD(1);
t278 = t239 * t254;
t277 = (qJD(5) * t238 + qJD(1)) * t234;
t276 = t259 * t238;
t264 = qJD(1) * t238;
t253 = t232 + t264;
t236 = sin(qJ(1));
t262 = qJD(3) * t236;
t257 = t235 * t262;
t274 = t253 * t239 - t257;
t272 = r_i_i_C(1) * t230;
t269 = t237 * pkin(5);
t267 = t232 * t235;
t246 = t254 * t236;
t224 = t230 * t246 - t274 * t231;
t225 = t274 * t230 + t231 * t246;
t266 = t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
t260 = qJD(3) * t239;
t256 = t235 * t260;
t243 = t253 * t236 + t256;
t226 = t230 * t278 + t243 * t231;
t227 = t243 * t230 - t231 * t278;
t265 = -t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
t263 = qJD(1) * t239;
t261 = qJD(3) * t238;
t251 = -qJD(5) - t264;
t250 = qJD(1) * t259;
t247 = -r_i_i_C(2) * t231 - t272;
t245 = -t247 + t255;
t244 = qJD(1) * t245;
t241 = qJD(2) - t249 * t238 + (-pkin(1) - pkin(7) - pkin(4) - t269) * qJD(1) + (t255 * t235 + t276) * qJD(3);
t228 = t261 * t271;
t1 = [t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t279 * t236 + t241 * t239, t263 (t239 * t250 + t245 * t262) * t238 + (t280 * t236 + t239 * t244) * t235, -t238 * t263 + t257 (t251 * t239 * t237 + (qJD(3) * t235 * t237 + t277) * t236) * pkin(5) + t266, t266; -t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t241 * t236 + t279 * t239, qJD(1) * t236 (t236 * t250 - t245 * t260) * t238 + (t236 * t244 - t280 * t239) * t235, -t236 * t264 - t256 (-t239 * t277 + (t251 * t236 - t256) * t237) * pkin(5) + t265, t265; 0, 0, t242 * t238 + (-t245 * t235 - t276) * qJD(3), t261, t228 + (t269 - t270) * t261 + (t247 * t232 - t258) * t235, -t267 * t272 + t228 + (-t230 * t261 - t231 * t267) * r_i_i_C(2);];
JaD_transl  = t1;
