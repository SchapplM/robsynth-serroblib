% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR14_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:12
% EndTime: 2019-02-26 22:38:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
t268 = sin(pkin(7));
t269 = sin(pkin(6));
t295 = t269 * t268;
t270 = cos(pkin(7));
t275 = cos(qJ(3));
t294 = t270 * t275;
t272 = sin(qJ(3));
t276 = cos(qJ(2));
t293 = t272 * t276;
t273 = sin(qJ(2));
t292 = t273 * t275;
t274 = sin(qJ(1));
t291 = t274 * t273;
t290 = t274 * t276;
t277 = cos(qJ(1));
t289 = t277 * t273;
t288 = t277 * t276;
t287 = qJD(1) * t269;
t286 = qJD(2) * t272;
t285 = t274 * t295;
t284 = t277 * t295;
t283 = t274 * t287;
t282 = t277 * t287;
t271 = cos(pkin(6));
t281 = t271 * t288 - t291;
t280 = -t271 * t290 - t289;
t279 = t271 * t289 + t290;
t278 = t271 * t291 - t288;
t267 = qJD(1) * t280 - qJD(2) * t279;
t266 = -qJD(1) * t281 + qJD(2) * t278;
t1 = [0, t282, -t266 * t268 + t270 * t282, -t266 * t294 + t280 * t286 + (-t272 * t279 - t275 * t284) * qJD(1) + (-t278 * t275 + (t270 * t280 + t285) * t272) * qJD(3), 0, 0; 0, t283, -t267 * t268 + t270 * t283, -t267 * t294 + t281 * t286 + (-t272 * t278 - t275 * t285) * qJD(1) + (t279 * t275 + (t270 * t281 - t284) * t272) * qJD(3), 0, 0; 0, 0, qJD(2) * t273 * t295, t271 * t268 * qJD(3) * t272 + ((t270 * t293 + t292) * qJD(3) + (t270 * t292 + t293) * qJD(2)) * t269, 0, 0;];
JgD_rot  = t1;
