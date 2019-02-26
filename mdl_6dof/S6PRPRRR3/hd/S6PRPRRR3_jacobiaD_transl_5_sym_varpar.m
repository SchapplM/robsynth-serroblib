% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:49
% EndTime: 2019-02-26 19:54:49
% DurationCPUTime: 0.11s
% Computational Cost: add. (256->46), mult. (358->85), div. (0->0), fcn. (331->11), ass. (0->43)
t258 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3);
t228 = pkin(12) + qJ(4);
t226 = qJ(5) + t228;
t222 = sin(t226);
t229 = qJD(4) + qJD(5);
t257 = t222 * t229;
t223 = cos(t226);
t256 = t223 * t229;
t230 = sin(pkin(11));
t231 = sin(pkin(6));
t255 = t230 * t231;
t232 = cos(pkin(11));
t254 = t231 * t232;
t234 = sin(qJ(2));
t253 = t231 * t234;
t233 = cos(pkin(6));
t252 = t233 * t234;
t235 = cos(qJ(2));
t251 = t233 * t235;
t216 = t230 * t235 + t232 * t252;
t246 = qJD(2) * t235;
t242 = t232 * t246;
t247 = qJD(2) * t234;
t243 = t230 * t247;
t211 = -t233 * t242 + t243;
t241 = t229 * t254 + t211;
t250 = (-t216 * t256 + t241 * t222) * r_i_i_C(1) + (t216 * t257 + t241 * t223) * r_i_i_C(2);
t218 = -t230 * t252 + t232 * t235;
t238 = t230 * t251 + t232 * t234;
t213 = t238 * qJD(2);
t240 = -t229 * t255 + t213;
t249 = (-t218 * t256 + t240 * t222) * r_i_i_C(1) + (t218 * t257 + t240 * t223) * r_i_i_C(2);
t244 = t231 * t246;
t237 = -t229 * t233 - t244;
t245 = t229 * t253;
t248 = (t237 * t222 - t223 * t245) * r_i_i_C(1) + (t222 * t245 + t237 * t223) * r_i_i_C(2);
t225 = cos(t228);
t239 = -r_i_i_C(1) * t223 + r_i_i_C(2) * t222 - pkin(4) * t225 - cos(pkin(12)) * pkin(3) - pkin(2);
t224 = sin(t228);
t236 = -pkin(4) * qJD(4) * t224 + (-r_i_i_C(1) * t222 - r_i_i_C(2) * t223) * t229;
t214 = -t233 * t243 + t242;
t212 = t216 * qJD(2);
t1 = [0, t218 * qJD(3) - t258 * t213 + t239 * t214 - t236 * t238, t214 (t213 * t224 + (-t218 * t225 - t224 * t255) * qJD(4)) * pkin(4) + t249, t249, 0; 0, t216 * qJD(3) - t258 * t211 + t236 * (-t230 * t234 + t232 * t251) + t239 * t212, t212 (t211 * t224 + (-t216 * t225 + t224 * t254) * qJD(4)) * pkin(4) + t250, t250, 0; 0 (qJD(3) * t234 + t236 * t235 + (t239 * t234 + t258 * t235) * qJD(2)) * t231, t231 * t247 (-t224 * t244 + (-t224 * t233 - t225 * t253) * qJD(4)) * pkin(4) + t248, t248, 0;];
JaD_transl  = t1;
