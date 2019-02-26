% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:12
% EndTime: 2019-02-26 20:48:13
% DurationCPUTime: 0.28s
% Computational Cost: add. (339->57), mult. (551->85), div. (0->0), fcn. (464->8), ass. (0->42)
t243 = sin(qJ(3));
t245 = cos(qJ(3));
t240 = pkin(9) + qJ(5);
t238 = sin(t240);
t239 = cos(t240);
t264 = qJD(6) * t238;
t275 = r_i_i_C(3) + qJ(6);
t277 = pkin(5) + r_i_i_C(1);
t249 = -t264 + (t277 * t238 - t275 * t239) * qJD(5);
t237 = cos(pkin(9)) * pkin(4) + pkin(3);
t251 = t275 * t238 + t277 * t239 + t237;
t276 = r_i_i_C(2) + pkin(8) + qJ(4);
t279 = t276 * t245;
t286 = (t243 * t251 - t279) * qJD(3) - qJD(4) * t243 + t249 * t245;
t258 = t276 * t243;
t283 = t251 * t245 + t258;
t282 = (-t237 * t243 - qJ(2) + t279) * qJD(1) + t239 * qJD(6);
t244 = sin(qJ(1));
t274 = t244 * t238;
t273 = t244 * t239;
t246 = cos(qJ(1));
t272 = t246 * t238;
t271 = qJD(1) * t244;
t270 = qJD(1) * t246;
t269 = qJD(3) * t243;
t268 = qJD(3) * t245;
t267 = qJD(3) * t246;
t265 = qJD(5) * t245;
t262 = t245 * qJD(4);
t261 = t246 * t243 * t239;
t260 = t245 * t267;
t255 = qJD(5) * t243 + qJD(1);
t254 = qJD(1) * t243 + qJD(5);
t253 = t243 * t273 + t272;
t250 = t244 * t268 + t254 * t246;
t248 = qJD(1) * t283;
t247 = t243 * t264 - t262 + qJD(2) + (t237 * t245 + t258) * qJD(3) + (-pkin(4) * sin(pkin(9)) - pkin(1) - pkin(7)) * qJD(1);
t232 = t250 * t239 - t255 * t274;
t231 = t250 * t238 + t255 * t273;
t230 = -t239 * t260 + (t243 * t272 + t273) * qJD(5) + t253 * qJD(1);
t229 = -qJD(5) * t261 - t238 * t260 - t239 * t270 + t254 * t274;
t1 = [-t275 * t229 - t277 * t230 + t282 * t244 + t247 * t246, t270, -t286 * t244 + t246 * t248, t244 * t269 - t245 * t270, t253 * qJD(6) - t277 * t231 + t275 * t232, t231; t275 * t231 + t277 * t232 + t247 * t244 - t282 * t246, t271, t244 * t248 + t286 * t246, -t243 * t267 - t245 * t271 -(t261 - t274) * qJD(6) + t275 * t230 - t277 * t229, t229; 0, 0, -t283 * qJD(3) + t249 * t243 + t262, t268 (-t275 * t265 + t277 * t269) * t238 + (-t275 * t269 + (-t277 * qJD(5) + qJD(6)) * t245) * t239, -t238 * t269 + t239 * t265;];
JaD_transl  = t1;
