% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:41
% EndTime: 2019-02-26 21:37:41
% DurationCPUTime: 0.22s
% Computational Cost: add. (322->48), mult. (303->61), div. (0->0), fcn. (221->10), ass. (0->37)
t234 = qJ(2) + pkin(10);
t231 = qJ(4) + t234;
t228 = cos(t231);
t267 = r_i_i_C(3) + qJ(5);
t277 = t267 * t228;
t222 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t234);
t214 = t222 * qJD(2);
t227 = sin(t231);
t226 = t227 * qJD(5);
t233 = qJD(2) + qJD(4);
t235 = sin(pkin(11));
t236 = cos(pkin(11));
t268 = r_i_i_C(2) * t235;
t275 = r_i_i_C(1) * t236 + pkin(4);
t245 = t275 - t268;
t276 = (-t245 * t227 + t277) * t233 + (t235 * r_i_i_C(1) + t236 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(3)) * qJD(1) + t214 + t226;
t260 = t233 * t268;
t274 = (qJD(5) + t260) * t228;
t238 = sin(qJ(1));
t265 = t228 * t238;
t240 = cos(qJ(1));
t264 = t233 * t240;
t263 = qJD(1) * t238;
t262 = qJD(1) * t240;
t258 = t227 * t263;
t257 = t227 * t262;
t255 = t267 * t227;
t253 = t267 * t238;
t250 = t274 * t240 + t275 * t258;
t249 = t275 * t233;
t248 = t275 * t240;
t247 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t234);
t246 = t274 * t238 + t257 * t268 + t262 * t277;
t243 = t226 + t233 * t277 + (-t249 + t260) * t227;
t242 = t247 * qJD(2) + (-t228 * t275 - t255) * t233;
t241 = qJD(3) + (-t245 * t228 - pkin(1) + t247 - t255) * qJD(1);
t1 = [-t276 * t238 + t241 * t240 (-t227 * t268 - t222 - t277) * t263 + t242 * t240 + t250, t262 (-t263 * t268 - t267 * t264) * t227 + (-qJD(1) * t253 - t233 * t248) * t228 + t250, t228 * t264 - t258, 0; t241 * t238 + t276 * t240 (-t227 * t275 + t222) * t262 + t242 * t238 + t246, t263, -t249 * t265 + (-qJD(1) * t248 - t233 * t253) * t227 + t246, t233 * t265 + t257, 0; 0, t214 + t243, 0, t243, t233 * t227, 0;];
JaD_transl  = t1;
