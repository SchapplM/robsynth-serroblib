% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRP11_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:42
% EndTime: 2019-02-26 21:13:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (71->39), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->42)
t262 = sin(pkin(7));
t266 = cos(pkin(6));
t291 = t262 * t266;
t263 = sin(pkin(6));
t269 = sin(qJ(1));
t290 = t263 * t269;
t272 = cos(qJ(1));
t289 = t263 * t272;
t265 = cos(pkin(7));
t268 = sin(qJ(3));
t288 = t265 * t268;
t261 = sin(pkin(12));
t287 = t269 * t261;
t264 = cos(pkin(12));
t286 = t269 * t264;
t285 = t272 * t261;
t284 = t272 * t264;
t283 = qJD(1) * t263;
t267 = sin(qJ(4));
t282 = qJD(3) * t267;
t281 = t269 * t283;
t280 = t272 * t283;
t279 = t262 * t281;
t278 = t262 * t280;
t257 = t266 * t284 - t287;
t277 = t257 * t265 - t262 * t289;
t259 = -t266 * t286 - t285;
t276 = t259 * t265 + t262 * t290;
t258 = t266 * t285 + t286;
t260 = -t266 * t287 + t284;
t271 = cos(qJ(3));
t275 = t258 * t271 + t277 * t268;
t274 = t260 * t271 + t276 * t268;
t273 = t268 * t291 + (t261 * t271 + t264 * t288) * t263;
t270 = cos(qJ(4));
t256 = t260 * qJD(1);
t255 = t259 * qJD(1);
t254 = t258 * qJD(1);
t253 = t257 * qJD(1);
t252 = -t255 * t262 + t265 * t281;
t251 = t253 * t262 + t265 * t280;
t1 = [0, 0, t251, -t254 * t268 + (t253 * t265 - t278) * t271 + t274 * qJD(3) (-t253 * t288 - t254 * t271 + t268 * t278) * t267 - t251 * t270 + (t274 * t270 + (-t259 * t262 + t265 * t290) * t267) * qJD(4) + (-t260 * t268 + t276 * t271) * t282, 0; 0, 0, t252, t256 * t268 + (-t255 * t265 - t279) * t271 + t275 * qJD(3) (t255 * t288 + t256 * t271 + t268 * t279) * t267 - t252 * t270 + (t275 * t270 + (-t257 * t262 - t265 * t289) * t267) * qJD(4) + (-t258 * t268 + t277 * t271) * t282, 0; 0, 0, 0, t273 * qJD(3) (t273 * t270 + (-t263 * t264 * t262 + t266 * t265) * t267) * qJD(4) + (t271 * t291 + (t264 * t265 * t271 - t261 * t268) * t263) * t282, 0;];
JgD_rot  = t1;
