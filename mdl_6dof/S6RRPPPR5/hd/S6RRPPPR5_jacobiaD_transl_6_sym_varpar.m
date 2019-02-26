% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:19
% EndTime: 2019-02-26 21:24:19
% DurationCPUTime: 0.39s
% Computational Cost: add. (234->72), mult. (756->115), div. (0->0), fcn. (683->8), ass. (0->49)
t241 = sin(qJ(2));
t238 = sin(pkin(9));
t239 = cos(pkin(9));
t240 = sin(qJ(6));
t243 = cos(qJ(6));
t274 = -pkin(5) - qJ(4);
t251 = t243 * r_i_i_C(1) - t240 * r_i_i_C(2) - t274;
t275 = -pkin(3) - qJ(5);
t252 = t240 * r_i_i_C(1) + t243 * r_i_i_C(2) - t275;
t277 = t251 * t238 + t252 * t239 + pkin(2);
t244 = cos(qJ(2));
t261 = r_i_i_C(3) + pkin(8) + pkin(4) + qJ(3);
t279 = t261 * t244;
t283 = t241 * t277 - t279;
t265 = t241 * qJD(3);
t282 = (-t241 * pkin(2) + t279) * qJD(2) + t265;
t254 = t238 * t240 - t239 * t243;
t255 = t238 * t243 + t239 * t240;
t278 = t254 * r_i_i_C(1) + t255 * r_i_i_C(2);
t281 = -t238 * qJD(4) - t239 * qJD(5) + qJD(6) * t278;
t242 = sin(qJ(1));
t273 = t242 * t244;
t245 = cos(qJ(1));
t272 = t245 * t238;
t271 = t245 * t239;
t270 = qJD(1) * t242;
t269 = qJD(1) * t245;
t268 = qJD(2) * t241;
t267 = qJD(2) * t244;
t266 = qJD(2) * t245;
t264 = t238 * t269;
t263 = t242 * t268;
t262 = t241 * t266;
t260 = t261 * t241;
t229 = t238 * t273 + t271;
t230 = t239 * t273 - t272;
t257 = t229 * t243 + t230 * t240;
t256 = t229 * t240 - t230 * t243;
t232 = t242 * t238 + t244 * t271;
t250 = -pkin(2) * t244 - pkin(1) - t260;
t246 = t244 * qJD(3) + t281 * t241 + (-t244 * t277 - t260) * qJD(2);
t231 = -t242 * t239 + t244 * t272;
t228 = t232 * qJD(1) - t239 * t263;
t227 = -t238 * t263 - t239 * t270 + t244 * t264;
t226 = -t264 + (t244 * t270 + t262) * t239;
t225 = t229 * qJD(1) + t238 * t262;
t224 = -t225 * t243 - t226 * t240 + (-t231 * t240 + t232 * t243) * qJD(6);
t223 = t225 * t240 - t226 * t243 + (-t231 * t243 - t232 * t240) * qJD(6);
t1 = [-t229 * qJD(4) - t230 * qJD(5) - t252 * t228 - t251 * t227 + (t256 * r_i_i_C(1) + t257 * r_i_i_C(2)) * qJD(6) - t282 * t242 + (-t242 * pkin(7) + t250 * t245) * qJD(1), t246 * t245 + t283 * t270, -t241 * t270 + t244 * t266, -t225, -t226, t223 * r_i_i_C(1) - t224 * r_i_i_C(2); t224 * r_i_i_C(1) + t223 * r_i_i_C(2) + t231 * qJD(4) + t232 * qJD(5) + t275 * t226 + t274 * t225 + t282 * t245 + (pkin(7) * t245 + t250 * t242) * qJD(1), t246 * t242 - t269 * t283, t241 * t269 + t242 * t267, t227, t228 (-t227 * t240 + t228 * t243) * r_i_i_C(1) + (-t227 * t243 - t228 * t240) * r_i_i_C(2) + (-t257 * r_i_i_C(1) + t256 * r_i_i_C(2)) * qJD(6); 0, -qJD(2) * t283 - t281 * t244 + t265, t268, t238 * t267, t239 * t267 (-t255 * r_i_i_C(1) + t254 * r_i_i_C(2)) * t241 * qJD(6) - t278 * t267;];
JaD_transl  = t1;
