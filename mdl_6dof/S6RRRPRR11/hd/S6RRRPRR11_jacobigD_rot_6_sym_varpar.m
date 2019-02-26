% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR11_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:51
% EndTime: 2019-02-26 22:21:51
% DurationCPUTime: 0.19s
% Computational Cost: add. (64->31), mult. (200->61), div. (0->0), fcn. (214->10), ass. (0->32)
t298 = qJD(5) - qJD(3);
t264 = sin(pkin(6));
t267 = sin(qJ(3));
t294 = t264 * t267;
t271 = cos(qJ(3));
t293 = t264 * t271;
t273 = cos(qJ(1));
t292 = t264 * t273;
t268 = sin(qJ(2));
t269 = sin(qJ(1));
t291 = t268 * t269;
t290 = t268 * t273;
t272 = cos(qJ(2));
t289 = t269 * t272;
t288 = t273 * t272;
t287 = qJD(1) * t264;
t286 = qJD(2) * t264;
t285 = t269 * t287;
t284 = t273 * t287;
t283 = t268 * t286;
t265 = cos(pkin(6));
t279 = t265 * t288 - t291;
t278 = t265 * t289 + t290;
t262 = t265 * t290 + t289;
t277 = t265 * t291 - t288;
t270 = cos(qJ(5));
t266 = sin(qJ(5));
t261 = -t277 * qJD(1) + t279 * qJD(2);
t260 = t278 * qJD(1) + t262 * qJD(2);
t259 = -t262 * qJD(1) - t278 * qJD(2);
t258 = -t279 * qJD(1) + t277 * qJD(2);
t1 = [0, t284, -t258, 0, t258 (t259 * t271 + t267 * t284) * t266 - (t259 * t267 - t271 * t284) * t270 - t298 * ((t267 * t277 + t269 * t293) * t266 - (t269 * t294 - t271 * t277) * t270); 0, t285, t260, 0, -t260 (t261 * t271 + t267 * t285) * t266 - (t261 * t267 - t271 * t285) * t270 + t298 * ((t262 * t267 + t271 * t292) * t266 + (t262 * t271 - t267 * t292) * t270); 0, 0, t283, 0, -t283 (t271 * t266 - t267 * t270) * t272 * t286 + t298 * ((-t265 * t271 + t268 * t294) * t266 + (t265 * t267 + t268 * t293) * t270);];
JgD_rot  = t1;
