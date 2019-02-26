% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:11
% EndTime: 2019-02-26 20:30:11
% DurationCPUTime: 0.25s
% Computational Cost: add. (341->52), mult. (380->78), div. (0->0), fcn. (298->9), ass. (0->44)
t226 = qJ(1) + pkin(9);
t222 = sin(t226);
t229 = sin(qJ(5));
t230 = cos(qJ(5));
t225 = pkin(10) + qJ(4);
t223 = cos(t225);
t250 = qJD(5) * t223;
t243 = -qJD(1) + t250;
t221 = sin(t225);
t253 = qJD(4) * t221;
t263 = -t229 * t253 + t243 * t230;
t269 = t263 * t222;
t258 = t230 * pkin(5);
t220 = pkin(4) + t258;
t249 = t221 * qJD(6);
t261 = pkin(5) * t229;
t257 = r_i_i_C(3) + qJ(6) + pkin(8);
t265 = t257 * t223;
t268 = (-t220 * t221 + t265) * qJD(4) - t250 * t261 + t249;
t260 = r_i_i_C(2) * t229;
t238 = r_i_i_C(1) * t230 + t220 - t260;
t232 = -t238 * t221 + t265;
t262 = pkin(5) + r_i_i_C(1);
t264 = r_i_i_C(2) * t230 + t262 * t229;
t255 = qJD(1) * t222;
t224 = cos(t226);
t254 = qJD(1) * t224;
t252 = qJD(4) * t223;
t251 = qJD(5) * t221;
t247 = t257 * t221;
t244 = pkin(7) + qJ(3) + t261;
t242 = qJD(1) * t223 - qJD(5);
t241 = qJD(5) * t258 + qJD(3);
t240 = t224 * t242;
t239 = t242 * t229;
t236 = t264 * t223;
t235 = -t220 * t223 - cos(pkin(10)) * pkin(3) - pkin(2) - t247;
t233 = t243 * t229 + t230 * t253;
t215 = t222 * t239 - t224 * t263;
t231 = qJD(6) * t223 + t264 * t251 + (-t238 * t223 - t247) * qJD(4);
t218 = t233 * t222 - t230 * t240;
t217 = t224 * t239 + t269;
t216 = t242 * t230 * t222 + t233 * t224;
t1 = [t218 * r_i_i_C(1) + t217 * r_i_i_C(2) + t241 * t224 - t268 * t222 + (-cos(qJ(1)) * pkin(1) - t244 * t222 + t235 * t224) * qJD(1), 0, t254, t231 * t224 - t232 * t255, t216 * r_i_i_C(2) + t262 * t215, -t221 * t255 + t224 * t252; -t216 * r_i_i_C(1) + t215 * r_i_i_C(2) + t241 * t222 + t268 * t224 + (-sin(qJ(1)) * pkin(1) + t244 * t224 + t235 * t222) * qJD(1), 0, t255, t231 * t222 + t232 * t254, -t217 * r_i_i_C(1) + t218 * r_i_i_C(2) + (-t229 * t240 - t269) * pkin(5), t221 * t254 + t222 * t252; 0, 0, 0, t232 * qJD(4) - qJD(5) * t236 + t249 (-t262 * t230 + t260) * t251 - qJD(4) * t236, t253;];
JaD_transl  = t1;
