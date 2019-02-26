% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:36
% EndTime: 2019-02-26 19:44:36
% DurationCPUTime: 0.20s
% Computational Cost: add. (168->46), mult. (442->88), div. (0->0), fcn. (445->11), ass. (0->40)
t259 = cos(pkin(6));
t254 = sin(pkin(11));
t257 = cos(pkin(11));
t261 = sin(qJ(2));
t262 = cos(qJ(2));
t265 = t262 * t254 + t261 * t257;
t240 = t265 * t259;
t272 = qJD(2) * t262;
t273 = qJD(2) * t261;
t279 = t254 * t273 - t257 * t272;
t278 = -r_i_i_C(3) - pkin(8) - qJ(4);
t277 = pkin(2) * qJD(2);
t255 = sin(pkin(10));
t256 = sin(pkin(6));
t276 = t255 * t256;
t258 = cos(pkin(10));
t275 = t256 * t258;
t274 = t259 * t261;
t253 = pkin(12) + qJ(5);
t251 = sin(t253);
t252 = cos(t253);
t269 = t251 * r_i_i_C(1) + t252 * r_i_i_C(2);
t235 = t279 * t259;
t242 = -t254 * t272 - t257 * t273;
t268 = t258 * t235 - t255 * t242;
t267 = t255 * t235 + t258 * t242;
t243 = t261 * t254 - t262 * t257;
t230 = t258 * t240 - t255 * t243;
t266 = t255 * t240 + t258 * t243;
t264 = t252 * r_i_i_C(1) - t251 * r_i_i_C(2) + cos(pkin(12)) * pkin(4) + pkin(3);
t263 = qJD(5) * t269;
t238 = t265 * t256;
t241 = t243 * qJD(2);
t239 = t243 * t259;
t236 = qJD(2) * t240;
t234 = qJD(2) * t238;
t233 = t279 * t256;
t226 = t255 * t236 + t258 * t241;
t223 = -t258 * t236 + t255 * t241;
t1 = [0, -t266 * qJD(4) - t278 * t267 - (t255 * t239 - t258 * t265) * t263 + (t255 * t274 - t258 * t262) * t277 + t264 * t226, 0, -t226, -t269 * t267 + ((-t251 * t276 + t252 * t266) * r_i_i_C(1) + (-t251 * t266 - t252 * t276) * r_i_i_C(2)) * qJD(5), 0; 0, t230 * qJD(4) + t278 * t268 - (-t258 * t239 - t255 * t265) * t263 + (-t255 * t262 - t258 * t274) * t277 + t264 * t223, 0, -t223, t269 * t268 + ((-t230 * t252 + t251 * t275) * r_i_i_C(1) + (t230 * t251 + t252 * t275) * r_i_i_C(2)) * qJD(5), 0; 0, t238 * qJD(4) + t278 * t233 - t264 * t234 + (-pkin(2) * t273 + t243 * t263) * t256, 0, t234, t269 * t233 + ((-t238 * t252 - t251 * t259) * r_i_i_C(1) + (t238 * t251 - t252 * t259) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
