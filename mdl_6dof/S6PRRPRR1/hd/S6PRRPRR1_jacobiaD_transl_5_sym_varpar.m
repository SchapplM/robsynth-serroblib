% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:06
% EndTime: 2019-02-26 20:04:06
% DurationCPUTime: 0.12s
% Computational Cost: add. (268->47), mult. (379->81), div. (0->0), fcn. (343->12), ass. (0->43)
t266 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8);
t236 = qJ(3) + pkin(12);
t233 = qJ(5) + t236;
t229 = sin(t233);
t235 = qJD(3) + qJD(5);
t265 = t229 * t235;
t230 = cos(t233);
t264 = t230 * t235;
t237 = sin(pkin(11));
t238 = sin(pkin(6));
t263 = t237 * t238;
t239 = cos(pkin(11));
t262 = t238 * t239;
t240 = cos(pkin(6));
t242 = sin(qJ(2));
t261 = t240 * t242;
t244 = cos(qJ(2));
t260 = t240 * t244;
t220 = t237 * t244 + t239 * t261;
t255 = qJD(2) * t244;
t252 = t239 * t255;
t256 = qJD(2) * t242;
t253 = t237 * t256;
t215 = -t240 * t252 + t253;
t250 = t235 * t262 + t215;
t259 = (-t220 * t264 + t250 * t229) * r_i_i_C(1) + (t220 * t265 + t250 * t230) * r_i_i_C(2);
t222 = -t237 * t261 + t239 * t244;
t247 = t237 * t260 + t239 * t242;
t217 = t247 * qJD(2);
t249 = -t235 * t263 + t217;
t258 = (-t222 * t264 + t249 * t229) * r_i_i_C(1) + (t222 * t265 + t249 * t230) * r_i_i_C(2);
t246 = -t235 * t240 - t238 * t255;
t254 = t235 * t238 * t242;
t257 = (t246 * t229 - t230 * t254) * r_i_i_C(1) + (t229 * t254 + t246 * t230) * r_i_i_C(2);
t226 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t236);
t251 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t236);
t248 = -r_i_i_C(1) * t230 + r_i_i_C(2) * t229 - pkin(2) + t251;
t223 = t226 * qJD(3);
t245 = t223 + (-r_i_i_C(1) * t229 - r_i_i_C(2) * t230) * t235;
t224 = t251 * qJD(3);
t218 = -t240 * t253 + t252;
t216 = t220 * qJD(2);
t1 = [0, t222 * qJD(4) - t266 * t217 + t248 * t218 - t245 * t247, -t217 * t226 + t222 * t224 + t223 * t263 + t258, t218, t258, 0; 0, t220 * qJD(4) - t266 * t215 + t245 * (-t237 * t242 + t239 * t260) + t248 * t216, -t215 * t226 + t220 * t224 - t223 * t262 + t259, t216, t259, 0; 0 (qJD(4) * t242 + t245 * t244 + (t248 * t242 + t266 * t244) * qJD(2)) * t238, t240 * t223 + (t224 * t242 + t226 * t255) * t238 + t257, t238 * t256, t257, 0;];
JaD_transl  = t1;
