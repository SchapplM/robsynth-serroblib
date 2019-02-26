% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR2
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
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:17
% EndTime: 2019-02-26 22:31:17
% DurationCPUTime: 0.24s
% Computational Cost: add. (492->57), mult. (398->69), div. (0->0), fcn. (286->10), ass. (0->49)
t238 = qJ(2) + qJ(3);
t235 = qJ(4) + t238;
t231 = cos(t235);
t279 = r_i_i_C(3) + qJ(5);
t291 = t279 * t231;
t230 = sin(t235);
t229 = t230 * qJD(5);
t236 = qJD(2) + qJD(3);
t232 = qJD(4) + t236;
t239 = sin(pkin(11));
t240 = cos(pkin(11));
t280 = r_i_i_C(2) * t239;
t289 = r_i_i_C(1) * t240 + pkin(4);
t253 = t289 - t280;
t241 = sin(qJ(2));
t278 = pkin(2) * qJD(2);
t268 = t241 * t278;
t233 = sin(t238);
t282 = pkin(3) * t236;
t271 = t233 * t282;
t290 = (-t253 * t230 + t291) * t232 + (r_i_i_C(1) * t239 + r_i_i_C(2) * t240 + pkin(7) + pkin(8) + pkin(9)) * qJD(1) - t268 - t271 + t229;
t269 = t232 * t280;
t288 = (qJD(5) + t269) * t231;
t234 = cos(t238);
t263 = t279 * t230;
t245 = (-t231 * t289 - t263) * t232 - t234 * t282;
t283 = pkin(3) * t233;
t242 = sin(qJ(1));
t276 = t231 * t242;
t244 = cos(qJ(1));
t275 = t232 * t244;
t274 = qJD(1) * t242;
t273 = qJD(1) * t244;
t266 = t230 * t273;
t264 = t230 * t274;
t261 = t279 * t242;
t258 = t288 * t244 + t289 * t264;
t257 = t289 * t230;
t256 = t289 * t232;
t255 = t289 * t244;
t254 = t288 * t242 + t266 * t280 + t273 * t291;
t251 = -t230 * t280 - t291;
t250 = t229 + t232 * t291 + (-t256 + t269) * t230;
t243 = cos(qJ(2));
t248 = -t243 * t278 + t245;
t247 = t250 - t271;
t246 = qJD(1) * (-pkin(2) * t243 - pkin(3) * t234 - t253 * t231 - pkin(1) - t263);
t226 = -pkin(2) * t241 - t283;
t1 = [-t290 * t242 + t244 * t246 (-t226 + t251) * t274 + t248 * t244 + t258 (t251 + t283) * t274 + t245 * t244 + t258 (-t274 * t280 - t279 * t275) * t230 + (-qJD(1) * t261 - t232 * t255) * t231 + t258, t231 * t275 - t264, 0; t242 * t246 + t290 * t244 (t226 - t257) * t273 + t248 * t242 + t254 (-t257 - t283) * t273 + t245 * t242 + t254, -t256 * t276 + (-qJD(1) * t255 - t232 * t261) * t230 + t254, t232 * t276 + t266, 0; 0, t247 - t268, t247, t250, t232 * t230, 0;];
JaD_transl  = t1;
