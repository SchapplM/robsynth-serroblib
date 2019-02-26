% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:42
% EndTime: 2019-02-26 20:43:42
% DurationCPUTime: 0.26s
% Computational Cost: add. (354->51), mult. (400->75), div. (0->0), fcn. (311->10), ass. (0->42)
t234 = cos(qJ(5));
t264 = pkin(5) * t234;
t222 = pkin(4) + t264;
t228 = qJ(3) + pkin(10);
t224 = sin(t228);
t226 = cos(t228);
t261 = r_i_i_C(3) + qJ(6) + pkin(8);
t243 = t261 * t226 - sin(qJ(3)) * pkin(3);
t232 = sin(qJ(5));
t250 = qJD(5) * t226 - qJD(1);
t247 = t250 * t232;
t253 = t224 * qJD(6);
t276 = (-t222 * t224 + t243) * qJD(3) - qJD(1) * (-qJ(4) - pkin(7)) + t253 - pkin(5) * t247;
t229 = qJ(1) + pkin(9);
t225 = sin(t229);
t256 = qJD(3) * t224;
t268 = -t232 * t256 + t250 * t234;
t275 = t268 * t225;
t263 = r_i_i_C(2) * t232;
t244 = r_i_i_C(1) * t234 + t222 - t263;
t238 = -t244 * t224 + t243;
t267 = pkin(5) + r_i_i_C(1);
t270 = -t261 * t224 - cos(qJ(3)) * pkin(3);
t269 = r_i_i_C(2) * t234 + t267 * t232;
t260 = pkin(1) * qJD(1);
t258 = qJD(1) * t225;
t227 = cos(t229);
t257 = qJD(1) * t227;
t255 = qJD(3) * t226;
t254 = qJD(5) * t224;
t249 = qJD(1) * t226 - qJD(5);
t248 = t227 * t249;
t246 = t249 * t232;
t241 = t226 * t269;
t239 = t234 * t256 + t247;
t237 = qJD(5) * t264 + qJD(4) + (-t222 * t226 - pkin(2) + t270) * qJD(1);
t218 = t225 * t246 - t227 * t268;
t236 = qJD(6) * t226 + t269 * t254 + (-t244 * t226 + t270) * qJD(3);
t221 = t239 * t225 - t234 * t248;
t220 = t227 * t246 + t275;
t219 = t249 * t234 * t225 + t239 * t227;
t1 = [-cos(qJ(1)) * t260 + t221 * r_i_i_C(1) + t220 * r_i_i_C(2) + t237 * t227 - t276 * t225, 0, t236 * t227 - t238 * t258, t257, t219 * r_i_i_C(2) + t267 * t218, -t224 * t258 + t227 * t255; -sin(qJ(1)) * t260 - t219 * r_i_i_C(1) + t218 * r_i_i_C(2) + t237 * t225 + t276 * t227, 0, t236 * t225 + t238 * t257, t258, -t220 * r_i_i_C(1) + t221 * r_i_i_C(2) + (-t232 * t248 - t275) * pkin(5), t224 * t257 + t225 * t255; 0, 0, t238 * qJD(3) - qJD(5) * t241 + t253, 0 (-t267 * t234 + t263) * t254 - qJD(3) * t241, t256;];
JaD_transl  = t1;
