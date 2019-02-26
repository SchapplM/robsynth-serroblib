% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:07
% EndTime: 2019-02-26 21:32:07
% DurationCPUTime: 0.26s
% Computational Cost: add. (217->60), mult. (636->95), div. (0->0), fcn. (594->8), ass. (0->43)
t240 = sin(qJ(5));
t243 = cos(qJ(5));
t270 = r_i_i_C(2) * t243;
t249 = r_i_i_C(1) * t240 + t270;
t271 = t249 * qJD(5) - qJD(3);
t269 = -pkin(4) - qJ(3);
t268 = pkin(8) - qJ(4);
t238 = sin(pkin(6));
t242 = sin(qJ(1));
t267 = t238 * t242;
t244 = cos(qJ(2));
t266 = t238 * t244;
t245 = cos(qJ(1));
t265 = t238 * t245;
t241 = sin(qJ(2));
t264 = t242 * t241;
t263 = t242 * t244;
t262 = t245 * t241;
t261 = t245 * t244;
t260 = qJD(1) * t242;
t259 = qJD(1) * t245;
t258 = qJD(2) * t241;
t257 = -r_i_i_C(3) - pkin(9) - pkin(3) - pkin(2);
t239 = cos(pkin(6));
t256 = t239 * t264;
t255 = t239 * t261;
t254 = t238 * t260;
t253 = t238 * t259;
t252 = t238 * t258;
t251 = qJD(2) * t239 + qJD(1);
t250 = -t243 * r_i_i_C(1) + t240 * r_i_i_C(2);
t230 = t239 * t263 + t262;
t248 = t239 * t262 + t263;
t247 = -t250 - t269;
t232 = t240 * t254;
t228 = -t255 + t264;
t227 = -qJD(1) * t256 - t242 * t258 + t251 * t261;
t226 = t230 * qJD(1) + t248 * qJD(2);
t225 = t248 * qJD(1) + t230 * qJD(2);
t224 = -qJD(1) * t255 - qJD(2) * t261 + t251 * t264;
t223 = -t240 * t253 - t224 * t243 + (-t230 * t240 - t243 * t267) * qJD(5);
t222 = -t243 * t253 + t224 * t240 + (-t230 * t243 + t240 * t267) * qJD(5);
t1 = [-pkin(1) * t259 + t232 * r_i_i_C(1) + t271 * t228 - t247 * t226 + ((t250 * qJD(5) - qJD(4)) * t245 + (-t268 + t270) * t260) * t238 + t257 * t227, -t271 * (-t256 + t261) - t247 * t225 - t257 * t224, -t224, -t253, t222 * r_i_i_C(1) - t223 * r_i_i_C(2), 0; -qJD(4) * t267 + t223 * r_i_i_C(1) + t222 * r_i_i_C(2) + t230 * qJD(3) + t269 * t224 + (-pkin(1) * t242 + t268 * t265) * qJD(1) + t257 * t225, t257 * t226 + t247 * t227 - t248 * t271, t226, -t254 (-t226 * t240 - t243 * t254) * r_i_i_C(1) + (-t226 * t243 + t232) * r_i_i_C(2) + ((-t228 * t243 - t240 * t265) * r_i_i_C(1) + (t228 * t240 - t243 * t265) * r_i_i_C(2)) * qJD(5), 0; 0 (-t271 * t241 + (t257 * t241 + t247 * t244) * qJD(2)) * t238, t252, 0, -t249 * t252 + ((t239 * t240 + t243 * t266) * r_i_i_C(1) + (t239 * t243 - t240 * t266) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
