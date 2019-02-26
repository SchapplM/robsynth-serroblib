% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobiaD_transl_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:08
% EndTime: 2019-02-26 21:37:08
% DurationCPUTime: 0.23s
% Computational Cost: add. (182->54), mult. (558->87), div. (0->0), fcn. (464->6), ass. (0->39)
t241 = sin(qJ(2));
t244 = cos(qJ(2));
t264 = qJD(2) * t244;
t243 = cos(qJ(4));
t261 = pkin(2) + pkin(8) + r_i_i_C(2);
t275 = t261 * qJD(2) + qJD(5) * t243 - qJD(3);
t277 = -qJ(3) * t264 + t275 * t241;
t242 = sin(qJ(1));
t258 = t242 * t264;
t245 = cos(qJ(1));
t266 = qJD(1) * t245;
t276 = t241 * t266 + t258;
t240 = sin(qJ(4));
t272 = r_i_i_C(3) + qJ(5);
t273 = -r_i_i_C(1) - pkin(4);
t249 = t273 * t240 + t272 * t243 - qJ(3);
t274 = pkin(3) + pkin(7);
t271 = t242 * t240;
t270 = t242 * t243;
t269 = t245 * t240;
t268 = t245 * t243;
t267 = qJD(1) * t242;
t265 = qJD(2) * t241;
t263 = qJD(4) * t244;
t262 = t240 * qJD(5);
t257 = t245 * t264;
t256 = qJD(4) * t241 + qJD(1);
t255 = qJD(1) * t241 + qJD(4);
t253 = t256 * t243;
t252 = t241 * t269 + t270;
t251 = -qJ(3) * t241 - t261 * t244 - pkin(1);
t248 = t249 * t244;
t247 = qJD(2) * t249;
t246 = (t272 * t240 - t273 * t243) * qJD(4) - t275;
t235 = t245 * t253 + (-t255 * t242 + t257) * t240;
t234 = -t243 * t257 + t252 * qJD(4) + (t241 * t270 + t269) * qJD(1);
t233 = t242 * t253 + (t255 * t245 + t258) * t240;
t232 = -qJD(4) * t268 - t276 * t243 + t256 * t271;
t1 = [t245 * t262 + t273 * t233 - t272 * t232 + t277 * t242 + (-t274 * t242 + t251 * t245) * qJD(1) (t245 * t247 + t261 * t267) * t241 + (t246 * t245 + t249 * t267) * t244, -t241 * t267 + t257, t252 * qJD(5) + t273 * t234 + t272 * t235, t234, 0; t242 * t262 - t273 * t235 + t272 * t234 - t277 * t245 + (t251 * t242 + t274 * t245) * qJD(1) (-t261 * t241 - t248) * t266 + (t241 * t247 + t246 * t244) * t242, t276 -(-t241 * t271 + t268) * qJD(5) + t272 * t233 + t273 * t232, t232, 0; 0, -qJD(2) * t248 + t246 * t241, t265 (-t272 * t263 - t273 * t265) * t243 + (t272 * t265 + (-t273 * qJD(4) - qJD(5)) * t244) * t240, -t240 * t263 - t243 * t265, 0;];
JaD_transl  = t1;
