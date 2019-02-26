% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:05
% EndTime: 2019-02-26 21:06:05
% DurationCPUTime: 0.26s
% Computational Cost: add. (168->55), mult. (518->88), div. (0->0), fcn. (433->6), ass. (0->39)
t241 = cos(qJ(3));
t237 = sin(qJ(4));
t240 = cos(qJ(4));
t267 = r_i_i_C(3) + qJ(5);
t269 = pkin(4) + r_i_i_C(1);
t245 = t267 * t237 + t269 * t240 + pkin(3);
t238 = sin(qJ(3));
t268 = pkin(8) + r_i_i_C(2);
t255 = t268 * t238;
t273 = t245 * t241 + t255;
t254 = t268 * t241;
t272 = -pkin(3) * t238 - qJ(2) + t254;
t257 = qJD(5) * t237;
t244 = -t257 + (t269 * t237 - t267 * t240) * qJD(4);
t270 = -pkin(1) - pkin(7);
t242 = cos(qJ(1));
t266 = t238 * t242;
t239 = sin(qJ(1));
t265 = t239 * t237;
t264 = t239 * t240;
t263 = qJD(1) * t239;
t262 = qJD(1) * t242;
t261 = qJD(3) * t238;
t260 = qJD(3) * t241;
t259 = qJD(3) * t242;
t258 = qJD(4) * t241;
t256 = qJD(5) * t240;
t253 = t240 * t266;
t252 = t241 * t259;
t250 = qJD(4) * t238 + qJD(1);
t249 = qJD(1) * t238 + qJD(4);
t248 = t249 * t242;
t247 = t237 * t242 + t238 * t264;
t243 = t238 * t257 + qJD(2) + (pkin(3) * t241 + t255) * qJD(3);
t232 = t240 * t248 + (-t250 * t237 + t240 * t260) * t239;
t231 = t250 * t264 + (t239 * t260 + t248) * t237;
t230 = -t240 * t252 + (t237 * t266 + t264) * qJD(4) + t247 * qJD(1);
t229 = -qJD(4) * t253 - t237 * t252 - t240 * t262 + t249 * t265;
t1 = [t239 * t256 - t269 * t230 - t267 * t229 + t243 * t242 + (t272 * t239 + t270 * t242) * qJD(1), t262, t273 * t262 + (-t244 * t241 + (-t238 * t245 + t254) * qJD(3)) * t239, t247 * qJD(5) - t269 * t231 + t267 * t232, t231, 0; -t242 * t256 + t269 * t232 + t267 * t231 + t243 * t239 + (t270 * t239 - t272 * t242) * qJD(1), t263 (t245 * t259 + t268 * t263) * t238 + (t245 * t263 + (-t268 * qJD(3) + t244) * t242) * t241 -(t253 - t265) * qJD(5) + t267 * t230 - t269 * t229, t229, 0; 0, 0, -t273 * qJD(3) + t244 * t238 (-t267 * t258 + t269 * t261) * t237 + (-t267 * t261 + (-t269 * qJD(4) + qJD(5)) * t241) * t240, -t237 * t261 + t240 * t258, 0;];
JaD_transl  = t1;
