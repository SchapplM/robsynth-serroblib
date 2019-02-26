% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:36
% EndTime: 2019-02-26 20:37:36
% DurationCPUTime: 0.23s
% Computational Cost: add. (270->57), mult. (436->85), div. (0->0), fcn. (344->8), ass. (0->49)
t240 = cos(qJ(5));
t231 = t240 * pkin(5) + pkin(4);
t238 = sin(qJ(4));
t241 = cos(qJ(4));
t237 = sin(qJ(5));
t278 = pkin(5) * t237;
t263 = qJD(5) * t278;
t273 = r_i_i_C(3) + pkin(9) + pkin(8);
t281 = t273 * t238;
t288 = (t231 * t241 + t281) * qJD(4) - t238 * t263 + qJD(3);
t236 = qJ(5) + qJ(6);
t233 = sin(t236);
t234 = cos(t236);
t287 = r_i_i_C(1) * t233 + r_i_i_C(2) * t234;
t286 = -r_i_i_C(1) * t234 + r_i_i_C(2) * t233;
t248 = t231 - t286;
t284 = t248 * t241 + t281;
t242 = cos(qJ(1));
t235 = qJD(5) + qJD(6);
t257 = t235 * t238 + qJD(1);
t283 = t242 * t257;
t280 = t287 * t235 + t263;
t268 = qJD(1) * t238;
t256 = t235 + t268;
t239 = sin(qJ(1));
t262 = qJD(4) * t239 * t241;
t279 = t256 * t242 + t262;
t265 = qJD(4) * t242;
t261 = t241 * t265;
t245 = t256 * t239 - t261;
t224 = t245 * t233 - t234 * t283;
t225 = t233 * t283 + t245 * t234;
t270 = t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
t251 = t257 * t239;
t226 = t279 * t233 + t234 * t251;
t227 = t233 * t251 - t279 * t234;
t269 = -t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
t267 = qJD(1) * t239;
t232 = qJD(1) * t242;
t266 = qJD(4) * t238;
t264 = qJD(5) * t240;
t258 = t273 * t241;
t255 = qJD(5) + t268;
t254 = pkin(7) - qJ(2) + t278;
t253 = -pkin(5) * t264 + qJD(2);
t250 = (-qJD(5) * t238 - qJD(1)) * t240;
t247 = t286 * t235 * t241 + t287 * t266;
t246 = -t231 * t238 - pkin(1) - qJ(3) + t258;
t1 = [t227 * r_i_i_C(1) + t226 * r_i_i_C(2) + t253 * t242 - t288 * t239 + (t254 * t239 + t246 * t242) * qJD(1), t232, -t267 (-t248 * t265 - t273 * t267) * t238 + (-t248 * t267 + (t273 * qJD(4) - t280) * t242) * t241 (t242 * t250 + (t255 * t239 - t261) * t237) * pkin(5) + t270, t270; -t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t253 * t239 + t288 * t242 + (t246 * t239 - t254 * t242) * qJD(1), t267, t232, t284 * t232 + (-t280 * t241 + (-t248 * t238 + t258) * qJD(4)) * t239 (t239 * t250 + (-t255 * t242 - t262) * t237) * pkin(5) + t269, t269; 0, 0, 0, -t284 * qJD(4) + t280 * t238 (t237 * t266 - t241 * t264) * pkin(5) + t247, t247;];
JaD_transl  = t1;
