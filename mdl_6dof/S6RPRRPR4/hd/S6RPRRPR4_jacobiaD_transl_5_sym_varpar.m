% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:48
% EndTime: 2019-02-26 21:02:48
% DurationCPUTime: 0.19s
% Computational Cost: add. (315->47), mult. (291->63), div. (0->0), fcn. (214->9), ass. (0->40)
t228 = pkin(10) + qJ(3);
t226 = qJ(4) + t228;
t223 = cos(t226);
t261 = r_i_i_C(3) + qJ(5);
t246 = t261 * t223;
t222 = sin(t226);
t221 = t222 * qJD(5);
t229 = qJD(3) + qJD(4);
t230 = sin(pkin(11));
t231 = cos(pkin(11));
t262 = r_i_i_C(2) * t230;
t269 = r_i_i_C(1) * t231 + pkin(4);
t239 = t269 - t262;
t224 = sin(t228);
t260 = pkin(3) * qJD(3);
t252 = t224 * t260;
t270 = (-t239 * t222 + t246) * t229 + (r_i_i_C(1) * t230 + r_i_i_C(2) * t231 + pkin(7) + pkin(8) + qJ(2)) * qJD(1) + t221 - t252;
t259 = t223 * t229;
t268 = qJD(5) * t223 + t259 * t262;
t264 = pkin(3) * t224;
t232 = sin(qJ(1));
t258 = t223 * t232;
t233 = cos(qJ(1));
t257 = t229 * t233;
t256 = qJD(1) * t232;
t255 = qJD(1) * t233;
t253 = t222 * t262;
t250 = t222 * t255;
t248 = t222 * t256;
t247 = t261 * t222;
t245 = t261 * t232;
t243 = t268 * t233 + t269 * t248;
t242 = t269 * t229;
t241 = t269 * t233;
t240 = t268 * t232 + t255 * t246 + t250 * t262;
t236 = -t222 * t242 + t229 * t253 + t261 * t259 + t221;
t225 = cos(t228);
t235 = -t225 * t260 + (-t223 * t269 - t247) * t229;
t234 = qJD(2) + (-t239 * t223 - pkin(3) * t225 - cos(pkin(10)) * pkin(2) - pkin(1) - t247) * qJD(1);
t1 = [-t270 * t232 + t234 * t233, t255 (-t246 - t253 + t264) * t256 + t235 * t233 + t243 (-t256 * t262 - t261 * t257) * t222 + (-qJD(1) * t245 - t229 * t241) * t223 + t243, t223 * t257 - t248, 0; t234 * t232 + t270 * t233, t256 (-t222 * t269 - t264) * t255 + t235 * t232 + t240, -t242 * t258 + (-qJD(1) * t241 - t229 * t245) * t222 + t240, t229 * t258 + t250, 0; 0, 0, t236 - t252, t236, t229 * t222, 0;];
JaD_transl  = t1;
