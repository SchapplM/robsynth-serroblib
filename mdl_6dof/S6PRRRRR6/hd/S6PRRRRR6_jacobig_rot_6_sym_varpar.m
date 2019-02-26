% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobig_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:52
% EndTime: 2019-02-26 20:21:53
% DurationCPUTime: 0.08s
% Computational Cost: add. (107->33), mult. (311->71), div. (0->0), fcn. (432->16), ass. (0->45)
t255 = sin(pkin(14));
t258 = sin(pkin(6));
t282 = t255 * t258;
t257 = sin(pkin(7));
t281 = t257 * t258;
t262 = cos(pkin(6));
t280 = t257 * t262;
t259 = cos(pkin(14));
t279 = t259 * t258;
t261 = cos(pkin(7));
t270 = cos(qJ(2));
t278 = t261 * t270;
t266 = sin(qJ(2));
t277 = t262 * t266;
t276 = t262 * t270;
t252 = t255 * t270 + t259 * t277;
t265 = sin(qJ(3));
t269 = cos(qJ(3));
t251 = -t255 * t266 + t259 * t276;
t272 = t251 * t261 - t257 * t279;
t242 = -t252 * t265 + t272 * t269;
t248 = -t251 * t257 - t261 * t279;
t256 = sin(pkin(8));
t260 = cos(pkin(8));
t275 = t242 * t260 + t248 * t256;
t254 = -t255 * t277 + t259 * t270;
t253 = -t255 * t276 - t259 * t266;
t271 = t253 * t261 + t255 * t281;
t244 = -t254 * t265 + t271 * t269;
t249 = -t253 * t257 + t261 * t282;
t274 = t244 * t260 + t249 * t256;
t246 = t269 * t280 + (-t265 * t266 + t269 * t278) * t258;
t250 = t262 * t261 - t270 * t281;
t273 = t246 * t260 + t250 * t256;
t268 = cos(qJ(4));
t267 = cos(qJ(5));
t264 = sin(qJ(4));
t263 = sin(qJ(5));
t247 = t265 * t280 + (t265 * t278 + t266 * t269) * t258;
t245 = t254 * t269 + t271 * t265;
t243 = t252 * t269 + t272 * t265;
t241 = -t246 * t256 + t250 * t260;
t240 = -t244 * t256 + t249 * t260;
t239 = -t242 * t256 + t248 * t260;
t1 = [0, t282, t249, t240, t245 * t264 - t274 * t268 (t245 * t268 + t274 * t264) * t263 - t240 * t267; 0, -t279, t248, t239, t243 * t264 - t275 * t268 (t243 * t268 + t275 * t264) * t263 - t239 * t267; 0, t262, t250, t241, t247 * t264 - t273 * t268 (t247 * t268 + t273 * t264) * t263 - t241 * t267;];
Jg_rot  = t1;
