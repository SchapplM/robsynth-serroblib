% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:18
% EndTime: 2019-02-26 20:24:18
% DurationCPUTime: 0.23s
% Computational Cost: add. (181->49), mult. (506->73), div. (0->0), fcn. (485->8), ass. (0->39)
t232 = sin(qJ(5));
t231 = sin(qJ(6));
t233 = cos(qJ(6));
t245 = t233 * r_i_i_C(1) - t231 * r_i_i_C(2);
t241 = pkin(5) + t245;
t234 = cos(qJ(5));
t259 = pkin(8) + r_i_i_C(3);
t248 = t259 * t234;
t236 = -t241 * t232 + t248;
t266 = t236 * qJD(5);
t249 = t259 * t232;
t237 = t241 * t234 + t249;
t255 = sin(pkin(9));
t256 = cos(pkin(9));
t257 = sin(qJ(1));
t258 = cos(qJ(1));
t226 = t258 * t255 - t257 * t256;
t263 = qJD(1) * t257;
t262 = qJD(1) * t258;
t261 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t260 = pkin(3) + pkin(7);
t254 = qJD(5) * t232;
t253 = qJD(5) * t234;
t252 = qJD(6) * t226;
t251 = qJD(6) * t232;
t250 = qJD(6) * t234;
t244 = t231 * r_i_i_C(1) + t233 * r_i_i_C(2);
t225 = -t257 * t255 - t258 * t256;
t223 = t225 * qJD(1);
t243 = t225 * t251 - t223;
t224 = t226 * qJD(1);
t242 = -t226 * t251 + t224;
t240 = qJD(6) * t244;
t239 = -qJD(6) * t225 + t223 * t232 + t226 * t253;
t238 = -t224 * t232 + t225 * t253 + t252;
t235 = t237 * qJD(5) - t232 * t240;
t222 = t242 * t231 + t239 * t233;
t221 = -t239 * t231 + t242 * t233;
t1 = [t245 * t252 - (-t244 - t260) * t223 - qJ(2) * t263 + (-qJ(4) + t236) * t224 + (qJD(4) + t235) * t225 + t261 * t258, t262, 0, t223, t237 * t223 + (-t234 * t240 + t266) * t226, t221 * r_i_i_C(1) - t222 * r_i_i_C(2); t222 * r_i_i_C(1) + t221 * r_i_i_C(2) + t260 * t224 + qJ(2) * t262 + (qJD(4) + (pkin(5) * t234 + t249) * qJD(5)) * t226 - (-pkin(5) * t232 - qJ(4) + t248) * t223 + t261 * t257, t263, 0, t224, t237 * t224 + (t244 * t250 - t266) * t225 (t243 * r_i_i_C(1) + t238 * r_i_i_C(2)) * t233 + (t238 * r_i_i_C(1) - t243 * r_i_i_C(2)) * t231; 0, 0, 0, 0, t235 (-t231 * t250 - t233 * t254) * r_i_i_C(2) + (-t231 * t254 + t233 * t250) * r_i_i_C(1);];
JaD_transl  = t1;
