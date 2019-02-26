% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:31
% EndTime: 2019-02-26 21:12:31
% DurationCPUTime: 0.21s
% Computational Cost: add. (265->53), mult. (428->81), div. (0->0), fcn. (338->8), ass. (0->47)
t235 = cos(qJ(4));
t227 = t235 * pkin(4) + pkin(3);
t233 = sin(qJ(3));
t236 = cos(qJ(3));
t266 = r_i_i_C(3) + pkin(9) + pkin(8);
t251 = t266 * t236;
t257 = qJD(4) * t235;
t280 = -pkin(4) * t257 + (-t227 * t233 - qJ(2) + t251) * qJD(1);
t231 = qJ(4) + qJ(5);
t228 = sin(t231);
t229 = cos(t231);
t276 = -r_i_i_C(1) * t229 + r_i_i_C(2) * t228;
t242 = t227 - t276;
t252 = t266 * t233;
t278 = t242 * t236 + t252;
t277 = r_i_i_C(1) * t228 + r_i_i_C(2) * t229;
t275 = t235 * (qJD(4) * t233 + qJD(1));
t230 = qJD(4) + qJD(5);
t232 = sin(qJ(4));
t271 = pkin(4) * t232;
t256 = qJD(4) * t271;
t240 = t277 * t230 + t256;
t234 = sin(qJ(1));
t262 = qJD(1) * t233;
t248 = t230 + t262;
t237 = cos(qJ(1));
t258 = qJD(3) * t237;
t253 = t236 * t258;
t273 = t248 * t234 - t253;
t254 = qJD(3) * t234 * t236;
t272 = t248 * t237 + t254;
t249 = -t230 * t233 - qJD(1);
t243 = t249 * t237;
t220 = t273 * t228 + t229 * t243;
t221 = t228 * t243 - t273 * t229;
t264 = -t220 * r_i_i_C(1) + t221 * r_i_i_C(2);
t244 = t249 * t234;
t222 = -t272 * t228 + t229 * t244;
t223 = t228 * t244 + t272 * t229;
t263 = t222 * r_i_i_C(1) - t223 * r_i_i_C(2);
t261 = qJD(1) * t234;
t260 = qJD(1) * t237;
t259 = qJD(3) * t233;
t246 = -qJD(4) - t262;
t241 = t276 * t230 * t236 + t277 * t259;
t239 = -t233 * t256 + qJD(2) + (t227 * t236 + t252) * qJD(3) + (-pkin(1) - pkin(7) - t271) * qJD(1);
t1 = [t221 * r_i_i_C(1) + t220 * r_i_i_C(2) + t280 * t234 + t239 * t237, t260, t278 * t260 + (-t240 * t236 + (-t242 * t233 + t251) * qJD(3)) * t234 (-t234 * t275 + (t246 * t237 - t254) * t232) * pkin(4) + t263, t263, 0; t223 * r_i_i_C(1) + t222 * r_i_i_C(2) + t239 * t234 - t280 * t237, t261 (t242 * t258 + t266 * t261) * t233 + (t242 * t261 + (-t266 * qJD(3) + t240) * t237) * t236 (t237 * t275 + (t246 * t234 + t253) * t232) * pkin(4) + t264, t264, 0; 0, 0, -t278 * qJD(3) + t240 * t233 (t232 * t259 - t236 * t257) * pkin(4) + t241, t241, 0;];
JaD_transl  = t1;
