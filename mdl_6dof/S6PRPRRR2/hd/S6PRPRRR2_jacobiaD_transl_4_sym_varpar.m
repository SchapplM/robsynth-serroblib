% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:22
% EndTime: 2019-02-26 19:54:22
% DurationCPUTime: 0.16s
% Computational Cost: add. (115->39), mult. (386->84), div. (0->0), fcn. (382->10), ass. (0->37)
t241 = sin(pkin(12));
t244 = cos(pkin(12));
t250 = cos(qJ(2));
t261 = qJD(2) * t250;
t248 = sin(qJ(2));
t262 = qJD(2) * t248;
t268 = t241 * t262 - t244 * t261;
t267 = -pkin(8) - r_i_i_C(3);
t266 = pkin(2) * qJD(2);
t243 = sin(pkin(6));
t247 = sin(qJ(4));
t265 = t243 * t247;
t249 = cos(qJ(4));
t264 = t243 * t249;
t246 = cos(pkin(6));
t263 = t246 * t248;
t258 = r_i_i_C(1) * t247 + r_i_i_C(2) * t249;
t227 = t268 * t246;
t234 = -t241 * t261 - t244 * t262;
t242 = sin(pkin(11));
t245 = cos(pkin(11));
t257 = t227 * t245 - t234 * t242;
t256 = t227 * t242 + t234 * t245;
t255 = t241 * t250 + t248 * t244;
t254 = t248 * t241 - t244 * t250;
t253 = r_i_i_C(1) * t249 - r_i_i_C(2) * t247 + pkin(3);
t252 = qJD(4) * t258;
t251 = qJD(2) * t255;
t233 = t254 * qJD(2);
t232 = t255 * t246;
t231 = t254 * t246;
t230 = t255 * t243;
t228 = t246 * t251;
t225 = t268 * t243;
t224 = -t232 * t242 - t245 * t254;
t222 = t232 * t245 - t242 * t254;
t1 = [0, -t267 * t256 - (t231 * t242 - t245 * t255) * t252 + (t242 * t263 - t245 * t250) * t266 + t253 * (t228 * t242 + t233 * t245) 0, -t258 * t256 + ((-t224 * t249 - t242 * t265) * r_i_i_C(1) + (t224 * t247 - t242 * t264) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t267 * t257 - (-t231 * t245 - t242 * t255) * t252 + (-t242 * t250 - t245 * t263) * t266 + t253 * (-t228 * t245 + t233 * t242) 0, t258 * t257 + ((-t222 * t249 + t245 * t265) * r_i_i_C(1) + (t222 * t247 + t245 * t264) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t267 * t225 + (-pkin(2) * t262 - t251 * t253 + t252 * t254) * t243, 0, t258 * t225 + ((-t230 * t249 - t246 * t247) * r_i_i_C(1) + (t230 * t247 - t246 * t249) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
