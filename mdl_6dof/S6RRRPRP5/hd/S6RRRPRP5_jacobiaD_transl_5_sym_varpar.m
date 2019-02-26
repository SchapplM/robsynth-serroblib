% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:32
% EndTime: 2019-02-26 22:11:32
% DurationCPUTime: 0.25s
% Computational Cost: add. (421->58), mult. (480->85), div. (0->0), fcn. (380->10), ass. (0->51)
t243 = qJ(3) + pkin(10);
t232 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t243);
t228 = t232 * qJD(3);
t233 = pkin(4) * cos(t243) + cos(qJ(3)) * pkin(3);
t231 = pkin(2) + t233;
t245 = sin(qJ(2));
t248 = cos(qJ(2));
t266 = t245 * qJD(4);
t277 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8);
t285 = t277 * t248;
t289 = (-t231 * t245 + t285) * qJD(2) + (pkin(7) + t232) * qJD(1) - t248 * t228 + t266;
t239 = qJ(5) + t243;
t235 = sin(t239);
t280 = r_i_i_C(2) * t235;
t236 = cos(t239);
t281 = r_i_i_C(1) * t236;
t255 = t231 - t280 + t281;
t251 = -t255 * t245 + t285;
t249 = cos(qJ(1));
t242 = qJD(3) + qJD(5);
t261 = t242 * t248 - qJD(1);
t287 = t249 * t261;
t279 = r_i_i_C(2) * t236;
t259 = r_i_i_C(1) * t235 + t279;
t284 = t259 * t242 + t228;
t246 = sin(qJ(1));
t271 = qJD(1) * t248;
t260 = -t242 + t271;
t269 = qJD(2) * t245;
t283 = -t246 * t269 + t260 * t249;
t275 = t242 * t245;
t267 = qJD(2) * t249;
t253 = t245 * t267 + t260 * t246;
t224 = t253 * t235 - t236 * t287;
t225 = t235 * t287 + t253 * t236;
t274 = t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
t257 = t261 * t246;
t226 = t283 * t235 + t236 * t257;
t227 = t235 * t257 - t283 * t236;
t273 = -t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
t272 = qJD(1) * t246;
t270 = qJD(1) * t249;
t268 = qJD(2) * t248;
t264 = t277 * t245;
t258 = t232 * t271 - t228;
t229 = t233 * qJD(3);
t254 = qJD(1) * t233 - t229 * t248 + t232 * t269;
t252 = t229 + (-t231 * t248 - pkin(1) - t264) * qJD(1);
t250 = qJD(4) * t248 + t284 * t245 + (-t255 * t248 - t264) * qJD(2);
t230 = t275 * t280;
t1 = [t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t289 * t246 + t252 * t249, t250 * t249 - t251 * t272, t258 * t246 + t254 * t249 + t274, -t245 * t272 + t248 * t267, t274, 0; -t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t252 * t246 + t289 * t249, t250 * t246 + t251 * t270, t254 * t246 - t258 * t249 + t273, t245 * t270 + t246 * t268, t273, 0; 0, t251 * qJD(2) - t284 * t248 + t266, t230 + (-t242 * t281 - t229) * t245 + (-t232 - t259) * t268, t269, -t268 * t279 + t230 + (-t235 * t268 - t236 * t275) * r_i_i_C(1), 0;];
JaD_transl  = t1;
