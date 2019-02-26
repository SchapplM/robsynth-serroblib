% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:22:27
% EndTime: 2019-02-26 21:22:27
% DurationCPUTime: 0.25s
% Computational Cost: add. (310->51), mult. (399->76), div. (0->0), fcn. (313->10), ass. (0->38)
t225 = qJ(2) + pkin(9);
t221 = sin(t225);
t223 = cos(t225);
t249 = pkin(3) + r_i_i_C(3) + pkin(8) + qJ(5);
t240 = t249 * t221 + sin(qJ(2)) * pkin(2);
t246 = pkin(5) * sin(pkin(10)) + qJ(4);
t250 = t223 * qJD(5);
t268 = (-t246 * t223 + t240) * qJD(2) - (cos(pkin(10)) * pkin(5) + pkin(4) + qJ(3) + pkin(7)) * qJD(1) - t221 * qJD(4) - t250;
t224 = pkin(10) + qJ(6);
t220 = sin(t224);
t222 = cos(t224);
t239 = r_i_i_C(1) * t220 + r_i_i_C(2) * t222 + t246;
t266 = -t239 * t223 + t240;
t230 = sin(qJ(1));
t245 = qJD(6) * t221 + qJD(1);
t265 = t230 * t245;
t232 = cos(qJ(1));
t264 = t232 * t245;
t261 = -t249 * t223 - cos(qJ(2)) * pkin(2);
t256 = qJD(1) * t230;
t255 = qJD(1) * t232;
t254 = qJD(2) * t221;
t253 = qJD(2) * t230;
t252 = qJD(2) * t232;
t251 = qJD(6) * t223;
t248 = t223 * t253;
t247 = t223 * t252;
t244 = -qJD(1) * t221 - qJD(6);
t238 = qJD(4) + (r_i_i_C(1) * t222 - r_i_i_C(2) * t220) * qJD(6);
t237 = t244 * t232 - t248;
t236 = t244 * t230 + t247;
t235 = qJD(3) + (-t246 * t221 - pkin(1) + t261) * qJD(1);
t233 = -qJD(5) * t221 + t238 * t223 + (-t239 * t221 + t261) * qJD(2);
t217 = t236 * t220 + t222 * t264;
t216 = -t220 * t264 + t236 * t222;
t215 = t237 * t220 - t222 * t265;
t214 = t220 * t265 + t237 * t222;
t1 = [t215 * r_i_i_C(1) + t214 * r_i_i_C(2) + t268 * t230 + t235 * t232, t233 * t232 + t266 * t256, t255, -t221 * t256 + t247, -t221 * t252 - t223 * t256, t216 * r_i_i_C(1) - t217 * r_i_i_C(2); t217 * r_i_i_C(1) + t216 * r_i_i_C(2) + t235 * t230 - t268 * t232, t233 * t230 - t255 * t266, t256, t221 * t255 + t248, -t221 * t253 + t223 * t255, -t214 * r_i_i_C(1) + t215 * r_i_i_C(2); 0, -qJD(2) * t266 + t238 * t221 + t250, 0, t254, qJD(2) * t223 (-t220 * t254 + t222 * t251) * r_i_i_C(2) + (t220 * t251 + t222 * t254) * r_i_i_C(1);];
JaD_transl  = t1;
