% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR5_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:11
% EndTime: 2019-02-26 20:21:11
% DurationCPUTime: 0.14s
% Computational Cost: add. (70->40), mult. (240->91), div. (0->0), fcn. (263->12), ass. (0->38)
t257 = sin(pkin(7));
t258 = sin(pkin(6));
t285 = t257 * t258;
t265 = cos(qJ(4));
t284 = t257 * t265;
t260 = cos(pkin(7));
t283 = t258 * t260;
t263 = sin(qJ(3));
t282 = t260 * t263;
t266 = cos(qJ(3));
t281 = t260 * t266;
t261 = cos(pkin(6));
t264 = sin(qJ(2));
t280 = t261 * t264;
t267 = cos(qJ(2));
t279 = t261 * t267;
t278 = t263 * t264;
t277 = t263 * t267;
t276 = t264 * t266;
t275 = t266 * t267;
t262 = sin(qJ(4));
t274 = qJD(3) * t262;
t256 = sin(pkin(13));
t259 = cos(pkin(13));
t252 = -t256 * t264 + t259 * t279;
t273 = t252 * t260 - t259 * t285;
t254 = -t256 * t279 - t259 * t264;
t272 = t254 * t260 + t256 * t285;
t253 = t256 * t267 + t259 * t280;
t271 = t256 * t280 - t259 * t267;
t270 = t260 * t277 + t276;
t269 = t253 * t266 + t273 * t263;
t268 = t272 * t263 - t266 * t271;
t251 = t271 * qJD(2);
t250 = t254 * qJD(2);
t249 = t253 * qJD(2);
t248 = t252 * qJD(2);
t1 = [0, 0, -t251 * t257, t268 * qJD(3) + t250 * t263 - t251 * t281 (t250 * t266 + t251 * t282) * t262 + t251 * t284 + (t268 * t265 + (-t254 * t257 + t256 * t283) * t262) * qJD(4) + (t263 * t271 + t272 * t266) * t274, 0; 0, 0, t249 * t257, t269 * qJD(3) + t248 * t263 + t249 * t281 (t248 * t266 - t249 * t282) * t262 - t249 * t284 + (t269 * t265 + (-t252 * t257 - t259 * t283) * t262) * qJD(4) + (-t253 * t263 + t273 * t266) * t274, 0; 0, 0, qJD(2) * t264 * t285, t261 * t257 * qJD(3) * t263 + (t270 * qJD(3) + (t260 * t276 + t277) * qJD(2)) * t258 (t257 * t266 * t274 + (t260 * t262 + t263 * t284) * qJD(4)) * t261 + ((-t257 * t267 * t262 + t270 * t265) * qJD(4) + (t260 * t275 - t278) * t274 + ((-t260 * t278 + t275) * t262 - t264 * t284) * qJD(2)) * t258, 0;];
JgD_rot  = t1;
