% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:31:56
% EndTime: 2019-02-26 20:31:56
% DurationCPUTime: 0.27s
% Computational Cost: add. (171->48), mult. (482->79), div. (0->0), fcn. (461->8), ass. (0->37)
t235 = sin(qJ(4));
t234 = sin(qJ(5));
t236 = cos(qJ(5));
t243 = r_i_i_C(1) * t236 - r_i_i_C(2) * t234 + pkin(4);
t237 = cos(qJ(4));
t260 = pkin(8) + r_i_i_C(3);
t250 = t260 * t237;
t239 = -t235 * t243 + t250;
t267 = t239 * qJD(4);
t249 = t260 * t235;
t266 = t243 * t237 + t249;
t256 = sin(pkin(9));
t257 = cos(pkin(9));
t258 = sin(qJ(1));
t259 = cos(qJ(1));
t229 = t256 * t259 - t257 * t258;
t263 = qJD(1) * t258;
t262 = qJD(1) * t259;
t261 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t255 = qJD(4) * t235;
t254 = qJD(4) * t237;
t253 = qJD(5) * t234;
t252 = qJD(5) * t236;
t251 = qJD(5) * t237;
t246 = r_i_i_C(1) * t234 + r_i_i_C(2) * t236;
t228 = -t256 * t258 - t257 * t259;
t226 = t228 * qJD(1);
t245 = t228 * t251 + t226;
t227 = t229 * qJD(1);
t244 = t229 * t251 + t227;
t242 = qJD(5) * t246;
t241 = qJD(5) * t228 + t226 * t237 - t229 * t255;
t240 = qJD(5) * t229 + t227 * t237 + t228 * t255;
t238 = qJD(4) * t266 - t235 * t242;
t225 = t234 * t245 + t236 * t240;
t224 = -t234 * t240 + t236 * t245;
t1 = [(-t227 * t234 + t228 * t252) * r_i_i_C(1) + (-t227 * t236 - t228 * t253) * r_i_i_C(2) - t227 * pkin(7) - qJ(2) * t263 - (-pkin(3) - t266) * t226 + (-t237 * t242 + t267) * t229 + t261 * t259, t262, 0, t227 * t239 + t228 * t238, r_i_i_C(1) * t224 - r_i_i_C(2) * t225, 0; t226 * pkin(7) + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + (pkin(4) * t235 - t250) * t228 * qJD(4) + qJ(2) * t262 + (pkin(4) * t237 + pkin(3) + t249) * t227 + t261 * t258, t263, 0, -t226 * t239 + t229 * t238 (r_i_i_C(1) * t244 + r_i_i_C(2) * t241) * t236 + (r_i_i_C(1) * t241 - r_i_i_C(2) * t244) * t234, 0; 0, 0, 0, t246 * t251 - t267 (-t235 * t253 + t236 * t254) * r_i_i_C(2) + (t234 * t254 + t235 * t252) * r_i_i_C(1), 0;];
JaD_transl  = t1;
