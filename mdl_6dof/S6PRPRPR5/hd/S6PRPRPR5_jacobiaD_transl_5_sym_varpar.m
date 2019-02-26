% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:46
% EndTime: 2019-02-26 19:48:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (215->42), mult. (429->77), div. (0->0), fcn. (411->9), ass. (0->38)
t262 = -pkin(4) + r_i_i_C(2);
t261 = r_i_i_C(1) + pkin(8) + qJ(3);
t260 = r_i_i_C(3) + qJ(5);
t237 = sin(pkin(10));
t238 = sin(pkin(6));
t259 = t237 * t238;
t239 = cos(pkin(10));
t258 = t238 * t239;
t242 = sin(qJ(2));
t257 = t238 * t242;
t240 = cos(pkin(6));
t256 = t240 * t242;
t243 = cos(qJ(2));
t255 = t240 * t243;
t254 = qJD(2) * t242;
t253 = qJD(2) * t243;
t252 = t238 * t253;
t251 = t237 * t254;
t250 = t239 * t253;
t228 = t237 * t243 + t239 * t256;
t236 = pkin(11) + qJ(4);
t234 = sin(t236);
t235 = cos(t236);
t249 = -t228 * t235 + t234 * t258;
t230 = -t237 * t256 + t239 * t243;
t248 = t230 * t235 + t234 * t259;
t247 = t240 * t234 + t235 * t257;
t246 = t237 * t255 + t239 * t242;
t245 = -t260 * t234 + t262 * t235 - cos(pkin(11)) * pkin(3) - pkin(2);
t244 = qJD(5) * t234 + (t262 * t234 + t260 * t235) * qJD(4);
t226 = -t240 * t251 + t250;
t225 = t246 * qJD(2);
t224 = t228 * qJD(2);
t223 = -t240 * t250 + t251;
t221 = t247 * qJD(4) + t234 * t252;
t219 = t248 * qJD(4) - t225 * t234;
t217 = -t249 * qJD(4) - t223 * t234;
t1 = [0, t230 * qJD(3) - t261 * t225 + t245 * t226 - t244 * t246, t226, t248 * qJD(5) + t260 * (-t225 * t235 + (-t230 * t234 + t235 * t259) * qJD(4)) + t262 * t219, t219, 0; 0, t228 * qJD(3) - t261 * t223 + t244 * (-t237 * t242 + t239 * t255) + t245 * t224, t224, -t249 * qJD(5) + t260 * (-t223 * t235 + (-t228 * t234 - t235 * t258) * qJD(4)) + t262 * t217, t217, 0; 0 (qJD(3) * t242 + t244 * t243 + (t245 * t242 + t261 * t243) * qJD(2)) * t238, t238 * t254, t247 * qJD(5) + t260 * (t235 * t252 + (-t234 * t257 + t235 * t240) * qJD(4)) + t262 * t221, t221, 0;];
JaD_transl  = t1;
