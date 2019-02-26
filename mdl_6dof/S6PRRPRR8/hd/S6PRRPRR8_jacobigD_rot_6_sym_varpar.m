% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR8_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:08
% EndTime: 2019-02-26 20:08:08
% DurationCPUTime: 0.21s
% Computational Cost: add. (70->40), mult. (240->92), div. (0->0), fcn. (263->12), ass. (0->40)
t263 = sin(qJ(3));
t266 = cos(qJ(3));
t256 = sin(pkin(12));
t259 = cos(pkin(12));
t267 = cos(qJ(2));
t261 = cos(pkin(6));
t264 = sin(qJ(2));
t278 = t261 * t264;
t269 = t256 * t278 - t259 * t267;
t277 = t261 * t267;
t254 = -t256 * t277 - t259 * t264;
t260 = cos(pkin(7));
t257 = sin(pkin(7));
t258 = sin(pkin(6));
t285 = t257 * t258;
t270 = t254 * t260 + t256 * t285;
t289 = t263 * t269 + t270 * t266;
t253 = t256 * t267 + t259 * t278;
t252 = -t256 * t264 + t259 * t277;
t271 = -t252 * t260 + t259 * t285;
t288 = t253 * t263 + t271 * t266;
t262 = sin(qJ(5));
t284 = t257 * t262;
t283 = t257 * t264;
t282 = t257 * t266;
t281 = t258 * t260;
t280 = t260 * t263;
t279 = t260 * t266;
t276 = t263 * t264;
t275 = t263 * t267;
t274 = t264 * t266;
t273 = t266 * t267;
t265 = cos(qJ(5));
t272 = qJD(3) * t265;
t268 = t260 * t273 - t276;
t251 = t269 * qJD(2);
t250 = t254 * qJD(2);
t249 = t253 * qJD(2);
t248 = t252 * qJD(2);
t1 = [0, 0, -t251 * t257, 0, t289 * qJD(3) + t250 * t266 + t251 * t280, -t251 * t284 - (t250 * t263 - t251 * t279) * t265 + ((-t254 * t257 + t256 * t281) * t265 - t289 * t262) * qJD(5) - (t270 * t263 - t266 * t269) * t272; 0, 0, t249 * t257, 0, -t288 * qJD(3) + t248 * t266 - t249 * t280, t249 * t284 - (t248 * t263 + t249 * t279) * t265 + ((-t252 * t257 - t259 * t281) * t265 + t288 * t262) * qJD(5) - (t253 * t266 - t271 * t263) * t272; 0, 0, t258 * qJD(2) * t283, 0, t261 * qJD(3) * t282 + (t268 * qJD(3) + (-t260 * t276 + t273) * qJD(2)) * t258 (-t257 * t263 * t272 + (t260 * t265 - t262 * t282) * qJD(5)) * t261 + ((-t257 * t267 * t265 - t268 * t262) * qJD(5) - (t260 * t275 + t274) * t272 + (t262 * t283 - (t260 * t274 + t275) * t265) * qJD(2)) * t258;];
JgD_rot  = t1;
