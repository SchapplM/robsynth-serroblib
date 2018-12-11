% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14_jacobig_rot_5_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobig_rot_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobig_rot_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:22
% EndTime: 2018-12-10 18:38:22
% DurationCPUTime: 0.11s
% Computational Cost: add. (194->45), mult. (194->65), div. (0->0), fcn. (204->25), ass. (0->50)
t279 = sin(pkin(6));
t286 = sin(qJ(1));
t290 = t286 * t279;
t288 = cos(qJ(1));
t289 = t288 * t279;
t287 = cos(qJ(2));
t285 = sin(qJ(2));
t284 = sin(qJ(4));
t283 = cos(pkin(6));
t282 = cos(pkin(7));
t281 = cos(pkin(8));
t280 = cos(pkin(14));
t278 = sin(pkin(7));
t277 = sin(pkin(8));
t276 = sin(pkin(14));
t275 = pkin(6) - qJ(2);
t274 = pkin(6) + qJ(2);
t273 = pkin(8) - qJ(4);
t272 = pkin(8) + qJ(4);
t271 = pkin(7) - pkin(14);
t270 = pkin(7) + pkin(14);
t269 = cos(t274);
t268 = sin(t275);
t267 = cos(t270);
t266 = sin(t271);
t265 = cos(t275) / 0.2e1;
t264 = sin(t274) / 0.2e1;
t263 = cos(t271) / 0.2e1;
t262 = sin(t270) / 0.2e1;
t261 = t265 - t269 / 0.2e1;
t260 = t265 + t269 / 0.2e1;
t259 = cos(t273) / 0.2e1 + cos(t272) / 0.2e1;
t258 = t264 - t268 / 0.2e1;
t257 = t264 + t268 / 0.2e1;
t256 = sin(t272) / 0.2e1 + sin(t273) / 0.2e1;
t255 = t263 - t267 / 0.2e1;
t254 = t263 + t267 / 0.2e1;
t253 = t262 - t266 / 0.2e1;
t252 = t262 + t266 / 0.2e1;
t251 = -t258 * t286 + t287 * t288;
t250 = -t260 * t286 - t285 * t288;
t249 = t258 * t288 + t286 * t287;
t248 = t260 * t288 - t285 * t286;
t247 = -t257 * t278 + t282 * t283;
t246 = -t250 * t278 + t282 * t290;
t245 = -t248 * t278 - t282 * t289;
t244 = t252 * t283 + t254 * t257 - t261 * t276;
t243 = t250 * t254 - t251 * t276 + t252 * t290;
t242 = t248 * t254 - t249 * t276 - t252 * t289;
t1 = [0, t290, 0, -t243 * t277 + t246 * t281 (t250 * t253 + t251 * t280 + t255 * t290) * t284 - t243 * t259 - t246 * t256, 0; 0, -t289, 0, -t242 * t277 + t245 * t281 (t248 * t253 + t249 * t280 - t255 * t289) * t284 - t242 * t259 - t245 * t256, 0; 1, t283, 0, -t244 * t277 + t247 * t281 (t253 * t257 + t255 * t283 + t261 * t280) * t284 - t244 * t259 - t247 * t256, 0;];
Jg_rot  = t1;
