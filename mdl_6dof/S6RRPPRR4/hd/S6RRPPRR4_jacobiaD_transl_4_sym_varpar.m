% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR4_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR4_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:12
% EndTime: 2019-02-26 21:30:12
% DurationCPUTime: 0.21s
% Computational Cost: add. (191->52), mult. (562->91), div. (0->0), fcn. (572->8), ass. (0->40)
t238 = cos(pkin(6));
t235 = sin(pkin(11));
t237 = cos(pkin(11));
t239 = sin(qJ(2));
t241 = cos(qJ(2));
t244 = t241 * t235 + t239 * t237;
t224 = t244 * t238;
t253 = qJD(2) * t241;
t251 = t238 * t253;
t254 = qJD(2) * t239;
t220 = t235 * t238 * t254 - t237 * t251;
t226 = -t235 * t253 - t237 * t254;
t240 = sin(qJ(1));
t242 = cos(qJ(1));
t228 = t239 * t235 - t241 * t237;
t245 = t240 * t224 + t242 * t228;
t264 = qJD(1) * t245 + t242 * t220 - t240 * t226;
t246 = t242 * t224 - t240 * t228;
t263 = qJD(1) * t246 - t240 * t220 - t242 * t226;
t262 = pkin(3) - r_i_i_C(2);
t261 = r_i_i_C(3) + qJ(4);
t260 = t239 * t240;
t259 = t239 * t242;
t258 = t240 * t241;
t257 = t241 * t242;
t256 = qJD(1) * t240;
t236 = sin(pkin(6));
t255 = qJD(2) * t236;
t252 = pkin(2) * t254;
t250 = -t238 * t239 * pkin(2) + (r_i_i_C(1) + pkin(8) + qJ(3)) * t236;
t223 = t228 * t238;
t247 = -t242 * t223 - t240 * t244;
t243 = t228 * qJD(2);
t234 = t241 * pkin(2) + pkin(1);
t227 = pkin(2) * t251 - t236 * qJD(3);
t221 = qJD(2) * t224;
t218 = t244 * t255;
t215 = t223 * t256 + (-qJD(1) * t244 - t221) * t242 + t240 * t243;
t212 = qJD(1) * t247 - t240 * t221 - t242 * t243;
t1 = [t247 * qJD(4) + t240 * t252 - t242 * t227 + t262 * t264 + t261 * t215 + (-t242 * t234 - t240 * t250) * qJD(1), -t245 * qJD(4) - t262 * t212 - t261 * t263 + ((t238 * t260 - t257) * qJD(2) + (-t238 * t257 + t260) * qJD(1)) * pkin(2), qJD(1) * t242 * t236, t212, 0, 0; -(t240 * t223 - t242 * t244) * qJD(4) - t242 * t252 - t240 * t227 - t262 * t263 + t261 * t212 + (-t240 * t234 + t242 * t250) * qJD(1), t246 * qJD(4) + t262 * t215 - t261 * t264 + ((-t238 * t259 - t258) * qJD(2) + (-t238 * t258 - t259) * qJD(1)) * pkin(2), t236 * t256, -t215, 0, 0; 0, -t261 * t228 * t255 - t262 * t218 + (qJD(4) * t244 - t252) * t236, 0, t218, 0, 0;];
JaD_transl  = t1;
