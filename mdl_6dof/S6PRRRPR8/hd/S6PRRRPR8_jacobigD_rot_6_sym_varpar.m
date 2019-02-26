% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR8_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:37
% EndTime: 2019-02-26 20:14:38
% DurationCPUTime: 0.16s
% Computational Cost: add. (70->40), mult. (240->92), div. (0->0), fcn. (263->12), ass. (0->40)
t259 = sin(pkin(7));
t260 = sin(pkin(6));
t289 = t259 * t260;
t264 = sin(qJ(4));
t288 = t259 * t264;
t265 = sin(qJ(3));
t287 = t259 * t265;
t266 = sin(qJ(2));
t286 = t259 * t266;
t262 = cos(pkin(7));
t285 = t260 * t262;
t284 = t262 * t265;
t268 = cos(qJ(3));
t283 = t262 * t268;
t263 = cos(pkin(6));
t282 = t263 * t266;
t269 = cos(qJ(2));
t281 = t263 * t269;
t280 = t265 * t266;
t279 = t265 * t269;
t278 = t266 * t268;
t277 = t268 * t269;
t267 = cos(qJ(4));
t276 = qJD(3) * t267;
t258 = sin(pkin(12));
t261 = cos(pkin(12));
t254 = -t258 * t266 + t261 * t281;
t275 = t254 * t262 - t261 * t289;
t256 = -t258 * t281 - t261 * t266;
t274 = t256 * t262 + t258 * t289;
t255 = t258 * t269 + t261 * t282;
t273 = t258 * t282 - t261 * t269;
t272 = t262 * t279 + t278;
t271 = t255 * t268 + t275 * t265;
t270 = t274 * t265 - t268 * t273;
t253 = t273 * qJD(2);
t252 = t256 * qJD(2);
t251 = t255 * qJD(2);
t250 = t254 * qJD(2);
t1 = [0, 0, -t253 * t259, t270 * qJD(3) + t252 * t265 - t253 * t283, 0 (t252 * t268 + t253 * t284) * t267 - t253 * t288 + (-t270 * t264 + (-t256 * t259 + t258 * t285) * t267) * qJD(4) + (t265 * t273 + t274 * t268) * t276; 0, 0, t251 * t259, t271 * qJD(3) + t250 * t265 + t251 * t283, 0 (t250 * t268 - t251 * t284) * t267 + t251 * t288 + (-t271 * t264 + (-t254 * t259 - t261 * t285) * t267) * qJD(4) + (-t255 * t265 + t275 * t268) * t276; 0, 0, t260 * qJD(2) * t286, t263 * qJD(3) * t287 + (t272 * qJD(3) + (t262 * t278 + t279) * qJD(2)) * t260, 0 (t259 * t268 * t276 + (t262 * t267 - t264 * t287) * qJD(4)) * t263 + ((-t259 * t269 * t267 - t272 * t264) * qJD(4) + (t262 * t277 - t280) * t276 + ((-t262 * t280 + t277) * t267 + t264 * t286) * qJD(2)) * t260;];
JgD_rot  = t1;
