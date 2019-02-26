% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:56
% EndTime: 2019-02-26 22:55:56
% DurationCPUTime: 0.20s
% Computational Cost: add. (145->52), mult. (459->96), div. (0->0), fcn. (432->10), ass. (0->34)
t237 = sin(pkin(14));
t238 = sin(pkin(7));
t240 = cos(pkin(14));
t265 = (r_i_i_C(1) * t237 + r_i_i_C(2) * t240) * t238 + pkin(10);
t264 = r_i_i_C(3) + qJ(3);
t242 = cos(pkin(6));
t245 = cos(qJ(2));
t246 = cos(qJ(1));
t253 = t246 * t245;
t243 = sin(qJ(2));
t244 = sin(qJ(1));
t256 = t244 * t243;
t247 = t242 * t256 - t253;
t250 = -t242 * t253 + t256;
t233 = t250 * qJD(1) + t247 * qJD(2);
t263 = t233 * t238;
t241 = cos(pkin(7));
t262 = t237 * t241;
t239 = sin(pkin(6));
t261 = t239 * t244;
t260 = t239 * t246;
t259 = t240 * t241;
t258 = t241 * t245;
t257 = t243 * t238;
t255 = t244 * t245;
t254 = t246 * t243;
t252 = qJD(1) * t239 * t241;
t249 = t242 * t255 + t254;
t248 = t242 * t254 + t255;
t236 = t247 * qJD(1) + t250 * qJD(2);
t235 = t249 * qJD(1) + t248 * qJD(2);
t234 = t248 * qJD(1) + t249 * qJD(2);
t232 = t246 * t252 - t263;
t1 = [(t235 * t262 + t236 * t240) * r_i_i_C(1) + (t235 * t259 - t236 * t237) * r_i_i_C(2) + t236 * pkin(2) + t241 * qJD(3) * t260 + (-t250 * qJD(3) - t264 * t235) * t238 + (-t246 * pkin(1) + (-t264 * t241 - t265) * t261) * qJD(1) (t233 * t240 + t234 * t262) * r_i_i_C(1) + (-t233 * t237 + t234 * t259) * r_i_i_C(2) + t233 * pkin(2) + (-t247 * qJD(3) - t264 * t234) * t238, t232, 0, 0, 0; (t233 * t262 - t234 * t240) * r_i_i_C(1) + (t233 * t259 + t234 * t237) * r_i_i_C(2) + t232 * r_i_i_C(3) - t234 * pkin(2) - qJ(3) * t263 + (t249 * t238 + t241 * t261) * qJD(3) + (-t244 * pkin(1) + (qJ(3) * t241 + t265) * t260) * qJD(1) (-t235 * t240 + t236 * t262) * r_i_i_C(1) + (t235 * t237 + t236 * t259) * r_i_i_C(2) - t235 * pkin(2) + (t248 * qJD(3) - t264 * t236) * t238, t235 * t238 + t244 * t252, 0, 0, 0; 0 (qJD(3) * t257 + ((-t237 * t258 - t240 * t243) * r_i_i_C(1) + (t237 * t243 - t240 * t258) * r_i_i_C(2) - t243 * pkin(2) + t264 * t245 * t238) * qJD(2)) * t239, t239 * qJD(2) * t257, 0, 0, 0;];
JaD_transl  = t1;
