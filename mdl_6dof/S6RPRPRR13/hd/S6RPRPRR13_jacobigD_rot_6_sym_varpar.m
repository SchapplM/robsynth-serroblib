% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR13_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:44
% EndTime: 2019-02-26 20:55:44
% DurationCPUTime: 0.18s
% Computational Cost: add. (71->39), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->42)
t257 = sin(pkin(12));
t259 = sin(pkin(6));
t260 = cos(pkin(12));
t264 = sin(qJ(3));
t267 = cos(qJ(3));
t261 = cos(pkin(7));
t283 = t261 * t267;
t258 = sin(pkin(7));
t262 = cos(pkin(6));
t286 = t258 * t262;
t291 = (-t257 * t264 + t260 * t283) * t259 + t267 * t286;
t268 = cos(qJ(1));
t279 = t268 * t260;
t265 = sin(qJ(1));
t282 = t265 * t257;
t256 = -t262 * t282 + t279;
t280 = t268 * t257;
t281 = t265 * t260;
t255 = -t262 * t281 - t280;
t285 = t259 * t265;
t270 = t255 * t261 + t258 * t285;
t290 = -t256 * t264 + t270 * t267;
t254 = t262 * t280 + t281;
t253 = t262 * t279 - t282;
t284 = t259 * t268;
t271 = -t253 * t261 + t258 * t284;
t289 = t254 * t264 + t271 * t267;
t278 = qJD(1) * t259;
t266 = cos(qJ(5));
t277 = qJD(3) * t266;
t275 = t265 * t278;
t274 = t268 * t278;
t273 = t258 * t275;
t272 = t258 * t274;
t263 = sin(qJ(5));
t252 = t256 * qJD(1);
t251 = t255 * qJD(1);
t250 = t254 * qJD(1);
t249 = t253 * qJD(1);
t248 = -t251 * t258 + t261 * t275;
t247 = t249 * t258 + t261 * t274;
t1 = [0, 0, t247, 0, -t250 * t267 + (-t249 * t261 + t272) * t264 + t290 * qJD(3), t247 * t263 - (t249 * t283 - t250 * t264 - t267 * t272) * t266 + ((-t255 * t258 + t261 * t285) * t266 - t290 * t263) * qJD(5) - (t256 * t267 + t270 * t264) * t277; 0, 0, t248, 0, t252 * t267 + (t251 * t261 + t273) * t264 - t289 * qJD(3), t248 * t263 - (-t251 * t283 + t252 * t264 - t267 * t273) * t266 + ((-t253 * t258 - t261 * t284) * t266 + t289 * t263) * qJD(5) - (t254 * t267 - t271 * t264) * t277; 0, 0, 0, 0, t291 * qJD(3) ((-t259 * t260 * t258 + t262 * t261) * t266 - t291 * t263) * qJD(5) - (t264 * t286 + (t260 * t261 * t264 + t257 * t267) * t259) * t277;];
JgD_rot  = t1;
