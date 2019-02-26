% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_jacobiaD_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:38
% EndTime: 2019-02-26 20:48:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (186->56), mult. (566->88), div. (0->0), fcn. (470->6), ass. (0->38)
t240 = sin(qJ(5));
t241 = sin(qJ(3));
t244 = cos(qJ(3));
t262 = pkin(3) + pkin(8) + r_i_i_C(2);
t278 = (-qJ(4) * t244 + t262 * t241 + qJ(2)) * qJD(1) + t240 * qJD(6);
t243 = cos(qJ(5));
t251 = t262 * qJD(3) + qJD(6) * t243 - qJD(4);
t273 = r_i_i_C(3) + qJ(6);
t274 = r_i_i_C(1) + pkin(5);
t247 = (t273 * t240 + t274 * t243) * qJD(5) - t251;
t242 = sin(qJ(1));
t245 = cos(qJ(1));
t267 = qJD(3) * t241;
t258 = t245 * t267;
t269 = qJD(1) * t244;
t277 = -t242 * t269 - t258;
t250 = t274 * t240 - t273 * t243 + qJ(4);
t272 = t242 * t244;
t271 = t245 * t240;
t270 = t245 * t243;
t268 = qJD(1) * t245;
t266 = qJD(3) * t242;
t265 = qJD(3) * t244;
t264 = qJD(5) * t241;
t260 = t244 * t268;
t259 = t241 * t266;
t256 = qJD(1) * t262;
t255 = qJD(5) * t244 + qJD(1);
t253 = t255 * t245;
t252 = t242 * t243 + t244 * t271;
t249 = qJD(1) * t250;
t248 = qJD(3) * t250;
t246 = qJ(4) * t267 + qJD(2) + (-pkin(1) - pkin(4) - pkin(7)) * qJD(1) + t251 * t244;
t233 = t243 * t253 + (-qJD(5) * t242 + t277) * t240;
t232 = t240 * t253 + (t258 + (qJD(5) + t269) * t242) * t243;
t231 = -t240 * t259 + (t243 * t272 + t271) * qJD(5) + t252 * qJD(1);
t230 = -t243 * t260 - qJD(5) * t270 + (t255 * t240 + t243 * t267) * t242;
t1 = [-t273 * t232 - t274 * t233 - t278 * t242 + t246 * t245, t268 (t245 * t256 + t250 * t266) * t244 + (t247 * t242 + t245 * t249) * t241, t259 - t260 -(t240 * t272 - t270) * qJD(6) - t273 * t231 + t274 * t230, -t230; -t273 * t230 - t274 * t231 + t246 * t242 + t278 * t245, qJD(1) * t242 (t242 * t256 - t245 * t248) * t244 + (t242 * t249 - t247 * t245) * t241, t277, t252 * qJD(6) - t274 * t232 + t273 * t233, t232; 0, 0, -t241 * t248 + t247 * t244, t265 (t273 * t264 + t274 * t265) * t243 + (t273 * t265 + (-t274 * qJD(5) + qJD(6)) * t241) * t240, t240 * t264 - t243 * t265;];
JaD_transl  = t1;
