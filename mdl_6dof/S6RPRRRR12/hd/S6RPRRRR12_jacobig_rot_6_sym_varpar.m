% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRRR12_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobig_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:10
% EndTime: 2019-02-26 21:21:10
% DurationCPUTime: 0.09s
% Computational Cost: add. (106->32), mult. (309->69), div. (0->0), fcn. (427->16), ass. (0->46)
t259 = sin(pkin(7));
t264 = cos(pkin(6));
t285 = t259 * t264;
t260 = sin(pkin(6));
t268 = sin(qJ(1));
t284 = t260 * t268;
t272 = cos(qJ(1));
t283 = t260 * t272;
t261 = cos(pkin(14));
t263 = cos(pkin(7));
t282 = t261 * t263;
t257 = sin(pkin(14));
t281 = t268 * t257;
t280 = t268 * t261;
t279 = t272 * t257;
t278 = t272 * t261;
t254 = t264 * t279 + t280;
t267 = sin(qJ(3));
t271 = cos(qJ(3));
t253 = t264 * t278 - t281;
t274 = t253 * t263 - t259 * t283;
t244 = -t254 * t267 + t274 * t271;
t250 = -t253 * t259 - t263 * t283;
t258 = sin(pkin(8));
t262 = cos(pkin(8));
t277 = t244 * t262 + t250 * t258;
t256 = -t264 * t281 + t278;
t255 = -t264 * t280 - t279;
t273 = t255 * t263 + t259 * t284;
t246 = -t256 * t267 + t273 * t271;
t251 = -t255 * t259 + t263 * t284;
t276 = t246 * t262 + t251 * t258;
t248 = t271 * t285 + (-t257 * t267 + t271 * t282) * t260;
t252 = -t260 * t261 * t259 + t264 * t263;
t275 = t248 * t262 + t252 * t258;
t270 = cos(qJ(4));
t269 = cos(qJ(5));
t266 = sin(qJ(4));
t265 = sin(qJ(5));
t249 = t267 * t285 + (t257 * t271 + t267 * t282) * t260;
t247 = t256 * t271 + t273 * t267;
t245 = t254 * t271 + t274 * t267;
t243 = -t248 * t258 + t252 * t262;
t242 = -t246 * t258 + t251 * t262;
t241 = -t244 * t258 + t250 * t262;
t1 = [0, 0, t251, t242, t247 * t266 - t276 * t270 (t247 * t270 + t276 * t266) * t265 - t242 * t269; 0, 0, t250, t241, t245 * t266 - t277 * t270 (t245 * t270 + t277 * t266) * t265 - t241 * t269; 1, 0, t252, t243, t249 * t266 - t275 * t270 (t249 * t270 + t275 * t266) * t265 - t243 * t269;];
Jg_rot  = t1;
