% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:51
% EndTime: 2019-02-26 19:45:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (132->49), mult. (436->94), div. (0->0), fcn. (430->10), ass. (0->40)
t277 = -pkin(2) - pkin(3);
t276 = pkin(8) + r_i_i_C(3);
t253 = sin(pkin(6));
t257 = sin(qJ(5));
t275 = t253 * t257;
t259 = cos(qJ(5));
t274 = t253 * t259;
t256 = cos(pkin(6));
t258 = sin(qJ(2));
t273 = t256 * t258;
t260 = cos(qJ(2));
t272 = t256 * t260;
t271 = qJD(2) * t258;
t270 = qJD(2) * t260;
t252 = sin(pkin(10));
t269 = t252 * t271;
t268 = t253 * t271;
t255 = cos(pkin(10));
t267 = t255 * t270;
t239 = -t256 * t267 + t269;
t244 = t252 * t260 + t255 * t273;
t240 = t244 * qJD(2);
t251 = sin(pkin(11));
t254 = cos(pkin(11));
t266 = -t239 * t254 + t240 * t251;
t245 = t252 * t272 + t255 * t258;
t241 = t245 * qJD(2);
t242 = -t256 * t269 + t267;
t265 = -t241 * t254 + t242 * t251;
t264 = t257 * r_i_i_C(1) + t259 * r_i_i_C(2);
t263 = t251 * t260 - t254 * t258;
t262 = t259 * r_i_i_C(1) - t257 * r_i_i_C(2) + pkin(4);
t261 = qJD(5) * t264;
t246 = -t252 * t273 + t255 * t260;
t243 = t252 * t258 - t255 * t272;
t238 = t263 * t253;
t235 = -t253 * t254 * t270 - t251 * t268;
t232 = t245 * t251 + t246 * t254;
t230 = t243 * t251 + t244 * t254;
t1 = [0, -t241 * qJ(3) + t246 * qJD(3) + t277 * t242 - t276 * t265 - (-t245 * t254 + t246 * t251) * t261 + t262 * (-t241 * t251 - t242 * t254) t242, 0, -t264 * t265 + ((-t232 * t259 + t252 * t275) * r_i_i_C(1) + (t232 * t257 + t252 * t274) * r_i_i_C(2)) * qJD(5), 0; 0, -t239 * qJ(3) + t244 * qJD(3) + t277 * t240 - t276 * t266 - (-t243 * t254 + t244 * t251) * t261 + t262 * (-t239 * t251 - t240 * t254) t240, 0, -t264 * t266 + ((-t230 * t259 - t255 * t275) * r_i_i_C(1) + (t230 * t257 - t255 * t274) * r_i_i_C(2)) * qJD(5), 0; 0, t276 * t235 + (-(t251 * t258 + t254 * t260) * t261 + qJD(3) * t258 + (qJ(3) * t260 + t277 * t258 + t262 * t263) * qJD(2)) * t253, t268, 0, t264 * t235 + ((t238 * t259 + t256 * t257) * r_i_i_C(1) + (-t238 * t257 + t256 * t259) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
