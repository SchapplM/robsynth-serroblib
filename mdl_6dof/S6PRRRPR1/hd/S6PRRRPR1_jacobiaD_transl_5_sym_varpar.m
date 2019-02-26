% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:39
% EndTime: 2019-02-26 20:10:39
% DurationCPUTime: 0.17s
% Computational Cost: add. (305->55), mult. (424->94), div. (0->0), fcn. (381->12), ass. (0->49)
t236 = qJ(3) + qJ(4);
t232 = sin(t236);
t268 = pkin(4) * t232;
t267 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8);
t266 = pkin(3) * qJD(3);
t231 = pkin(12) + t236;
t229 = sin(t231);
t235 = qJD(3) + qJD(4);
t265 = t229 * t235;
t230 = cos(t231);
t264 = t230 * t235;
t233 = cos(t236);
t263 = t233 * t235;
t237 = sin(pkin(11));
t238 = sin(pkin(6));
t262 = t237 * t238;
t239 = cos(pkin(11));
t261 = t238 * t239;
t240 = cos(pkin(6));
t242 = sin(qJ(2));
t260 = t240 * t242;
t244 = cos(qJ(2));
t259 = t240 * t244;
t220 = t237 * t244 + t239 * t260;
t254 = qJD(2) * t244;
t251 = t239 * t254;
t255 = qJD(2) * t242;
t252 = t237 * t255;
t215 = -t240 * t251 + t252;
t250 = t235 * t261 + t215;
t258 = (-t220 * t264 + t250 * t229) * r_i_i_C(1) + (t220 * t265 + t250 * t230) * r_i_i_C(2);
t222 = -t237 * t260 + t239 * t244;
t247 = t237 * t259 + t239 * t242;
t217 = t247 * qJD(2);
t249 = -t235 * t262 + t217;
t257 = (-t222 * t264 + t249 * t229) * r_i_i_C(1) + (t222 * t265 + t249 * t230) * r_i_i_C(2);
t246 = -t235 * t240 - t238 * t254;
t253 = t235 * t238 * t242;
t256 = (t246 * t229 - t230 * t253) * r_i_i_C(1) + (t229 * t253 + t246 * t230) * r_i_i_C(2);
t243 = cos(qJ(3));
t248 = -t243 * pkin(3) - pkin(4) * t233 - r_i_i_C(1) * t230 + r_i_i_C(2) * t229 - pkin(2);
t241 = sin(qJ(3));
t223 = -t235 * t268 - t241 * t266;
t245 = t223 + (-r_i_i_C(1) * t229 - r_i_i_C(2) * t230) * t235;
t226 = -t241 * pkin(3) - t268;
t224 = -pkin(4) * t263 - t243 * t266;
t218 = -t240 * t252 + t251;
t216 = t220 * qJD(2);
t1 = [0, t222 * qJD(5) - t267 * t217 + t248 * t218 - t245 * t247, -t217 * t226 + t222 * t224 + t223 * t262 + t257 (-t222 * t263 + t249 * t232) * pkin(4) + t257, t218, 0; 0, t220 * qJD(5) - t267 * t215 + t245 * (-t237 * t242 + t239 * t259) + t248 * t216, -t215 * t226 + t220 * t224 - t223 * t261 + t258 (-t220 * t263 + t250 * t232) * pkin(4) + t258, t216, 0; 0 (qJD(5) * t242 + t245 * t244 + (t248 * t242 + t267 * t244) * qJD(2)) * t238, t240 * t223 + (t224 * t242 + t226 * t254) * t238 + t256 (t246 * t232 - t233 * t253) * pkin(4) + t256, t238 * t255, 0;];
JaD_transl  = t1;
