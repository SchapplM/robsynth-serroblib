% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:13
% EndTime: 2019-02-26 19:45:13
% DurationCPUTime: 0.16s
% Computational Cost: add. (142->40), mult. (472->79), div. (0->0), fcn. (472->10), ass. (0->33)
t243 = sin(qJ(5));
t245 = cos(qJ(5));
t253 = t245 * r_i_i_C(1) - t243 * r_i_i_C(2);
t260 = t253 * qJD(5) + qJD(4);
t259 = pkin(2) * qJD(2);
t239 = sin(pkin(6));
t258 = t239 * t243;
t257 = t239 * t245;
t242 = cos(pkin(6));
t244 = sin(qJ(2));
t256 = t242 * t244;
t254 = pkin(3) + pkin(8) + r_i_i_C(3);
t237 = sin(pkin(11));
t240 = cos(pkin(11));
t246 = cos(qJ(2));
t252 = t246 * t237 + t244 * t240;
t235 = t244 * t237 - t246 * t240;
t251 = t243 * r_i_i_C(1) + t245 * r_i_i_C(2) + qJ(4);
t250 = t252 * t242;
t234 = t252 * qJD(2);
t233 = t235 * qJD(2);
t249 = qJD(2) * t250;
t248 = t242 * t233;
t241 = cos(pkin(10));
t238 = sin(pkin(10));
t232 = t235 * t242;
t230 = t235 * t239;
t228 = t239 * t234;
t224 = -t238 * t232 + t241 * t252;
t222 = t241 * t232 + t238 * t252;
t219 = t241 * t233 + t238 * t249;
t217 = t238 * t233 - t241 * t249;
t1 = [0 (t238 * t256 - t241 * t246) * t259 - t260 * (t241 * t235 + t238 * t250) - t251 * (t241 * t234 - t238 * t248) + t254 * t219, 0, -t219, -t253 * t219 + ((-t224 * t243 - t238 * t257) * r_i_i_C(1) + (-t224 * t245 + t238 * t258) * r_i_i_C(2)) * qJD(5), 0; 0 (-t238 * t246 - t241 * t256) * t259 - t260 * (t238 * t235 - t241 * t250) - t251 * (t238 * t234 + t241 * t248) + t254 * t217, 0, -t217, -t253 * t217 + ((-t222 * t243 + t241 * t257) * r_i_i_C(1) + (-t222 * t245 - t241 * t258) * r_i_i_C(2)) * qJD(5), 0; 0, -t254 * t228 + (t260 * t252 + (-t244 * pkin(2) - t251 * t235) * qJD(2)) * t239, 0, t228, t253 * t228 + ((-t230 * t243 - t242 * t245) * r_i_i_C(1) + (-t230 * t245 + t242 * t243) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
