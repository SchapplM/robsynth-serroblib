% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:47
% EndTime: 2019-02-26 19:58:47
% DurationCPUTime: 0.18s
% Computational Cost: add. (233->52), mult. (490->93), div. (0->0), fcn. (464->10), ass. (0->41)
t268 = -pkin(4) + r_i_i_C(2);
t267 = r_i_i_C(1) + qJ(4) + pkin(8);
t266 = r_i_i_C(3) + qJ(5);
t240 = sin(pkin(10));
t241 = sin(pkin(6));
t265 = t240 * t241;
t242 = cos(pkin(10));
t264 = t241 * t242;
t245 = sin(qJ(3));
t263 = t241 * t245;
t246 = sin(qJ(2));
t262 = t241 * t246;
t243 = cos(pkin(6));
t261 = t243 * t246;
t248 = cos(qJ(2));
t260 = t243 * t248;
t259 = qJD(2) * t246;
t258 = qJD(2) * t248;
t257 = t241 * t258;
t256 = t240 * t259;
t255 = t242 * t258;
t231 = t240 * t248 + t242 * t261;
t239 = qJ(3) + pkin(11);
t237 = sin(t239);
t238 = cos(t239);
t254 = -t231 * t238 + t237 * t264;
t233 = -t240 * t261 + t242 * t248;
t253 = t233 * t238 + t237 * t265;
t252 = t243 * t237 + t238 * t262;
t251 = t240 * t260 + t242 * t246;
t247 = cos(qJ(3));
t250 = -t247 * pkin(3) - t266 * t237 + t268 * t238 - pkin(2);
t249 = t237 * qJD(5) + (-pkin(3) * t245 + t268 * t237 + t266 * t238) * qJD(3);
t229 = -t243 * t256 + t255;
t228 = t251 * qJD(2);
t227 = t231 * qJD(2);
t226 = -t243 * t255 + t256;
t224 = qJD(3) * t252 + t237 * t257;
t222 = t253 * qJD(3) - t228 * t237;
t220 = -t254 * qJD(3) - t226 * t237;
t1 = [0, t233 * qJD(4) - t267 * t228 + t250 * t229 - t249 * t251, t253 * qJD(5) + t266 * (-t228 * t238 + (-t233 * t237 + t238 * t265) * qJD(3)) + t268 * t222 + (t228 * t245 + (-t233 * t247 - t240 * t263) * qJD(3)) * pkin(3), t229, t222, 0; 0, t231 * qJD(4) - t267 * t226 + t250 * t227 + t249 * (-t240 * t246 + t242 * t260) -t254 * qJD(5) + t266 * (-t226 * t238 + (-t231 * t237 - t238 * t264) * qJD(3)) + t268 * t220 + (t226 * t245 + (-t231 * t247 + t242 * t263) * qJD(3)) * pkin(3), t227, t220, 0; 0 ((qJD(2) * t250 + qJD(4)) * t246 + (t267 * qJD(2) + t249) * t248) * t241, t252 * qJD(5) + t266 * (t238 * t257 + (-t237 * t262 + t238 * t243) * qJD(3)) + t268 * t224 + (-t245 * t257 + (-t243 * t245 - t247 * t262) * qJD(3)) * pkin(3), t241 * t259, t224, 0;];
JaD_transl  = t1;
