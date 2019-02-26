% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:02
% EndTime: 2019-02-26 21:05:03
% DurationCPUTime: 0.26s
% Computational Cost: add. (425->58), mult. (488->84), div. (0->0), fcn. (386->10), ass. (0->51)
t244 = qJ(4) + pkin(10);
t234 = pkin(5) * cos(t244) + cos(qJ(4)) * pkin(4);
t228 = t234 * qJD(4);
t230 = pkin(3) + t234;
t246 = sin(qJ(3));
t249 = cos(qJ(3));
t280 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
t287 = t280 * t249;
t297 = (-t230 * t246 - qJ(2) + t287) * qJD(1) - t228;
t233 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t244);
t227 = t233 * qJD(4);
t243 = qJD(4) + qJD(6);
t240 = qJ(6) + t244;
t236 = sin(t240);
t237 = cos(t240);
t294 = r_i_i_C(1) * t236 + r_i_i_C(2) * t237;
t254 = t294 * t243 + t227;
t282 = r_i_i_C(2) * t236;
t283 = r_i_i_C(1) * t237;
t255 = t230 - t282 + t283;
t296 = (-t255 * t246 + t287) * qJD(3) + qJD(5) * t246 - t254 * t249;
t264 = t280 * t246;
t295 = t255 * t249 + t264;
t247 = sin(qJ(1));
t277 = qJD(1) * t246;
t260 = t243 + t277;
t250 = cos(qJ(1));
t272 = qJD(3) * t250;
t286 = t260 * t247 - t249 * t272;
t273 = qJD(3) * t249;
t285 = t247 * t273 + t260 * t250;
t261 = -t243 * t246 - qJD(1);
t256 = t261 * t250;
t223 = t286 * t236 + t237 * t256;
t224 = t236 * t256 - t286 * t237;
t279 = -t223 * r_i_i_C(1) + t224 * r_i_i_C(2);
t257 = t261 * t247;
t225 = -t285 * t236 + t237 * t257;
t226 = t236 * t257 + t285 * t237;
t278 = t225 * r_i_i_C(1) - t226 * r_i_i_C(2);
t276 = qJD(1) * t247;
t275 = qJD(1) * t250;
t274 = qJD(3) * t246;
t270 = t249 * qJD(5);
t269 = t243 * t283;
t268 = t249 * t243 * t282 + t294 * t274;
t258 = -t233 * t277 - t227;
t253 = qJD(1) * t234 + t228 * t246 + t233 * t273;
t252 = qJD(1) * t295;
t251 = -t270 - t246 * t227 + qJD(2) + (t230 * t249 + t264) * qJD(3) + (-pkin(1) - pkin(7) - t233) * qJD(1);
t1 = [t224 * r_i_i_C(1) + t223 * r_i_i_C(2) + t297 * t247 + t251 * t250, t275, t296 * t247 + t250 * t252, -t253 * t247 + t258 * t250 + t278, t247 * t274 - t249 * t275, t278; t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t251 * t247 - t297 * t250, t276, t247 * t252 - t296 * t250, t258 * t247 + t253 * t250 + t279, -t246 * t272 - t249 * t276, t279; 0, 0, -t295 * qJD(3) + t254 * t246 + t270, t233 * t274 + (-t228 - t269) * t249 + t268, t273, -t249 * t269 + t268;];
JaD_transl  = t1;
