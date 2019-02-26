% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:43
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.13s
% Computational Cost: add. (66->25), mult. (181->56), div. (0->0), fcn. (208->12), ass. (0->34)
t252 = sin(pkin(12));
t258 = cos(pkin(6));
t272 = t252 * t258;
t253 = sin(pkin(7));
t254 = sin(pkin(6));
t271 = t253 * t254;
t270 = t253 * t258;
t257 = cos(pkin(7));
t269 = t254 * t257;
t255 = cos(pkin(13));
t268 = t255 * t257;
t256 = cos(pkin(12));
t267 = t256 * t258;
t250 = qJ(4) + qJ(5);
t247 = sin(t250);
t266 = qJD(3) * t247;
t251 = sin(pkin(13));
t243 = -t251 * t252 + t255 * t267;
t265 = t243 * t257 - t256 * t271;
t245 = -t251 * t256 - t255 * t272;
t264 = t245 * t257 + t252 * t271;
t244 = t251 * t267 + t252 * t255;
t259 = sin(qJ(3));
t260 = cos(qJ(3));
t263 = t244 * t260 + t259 * t265;
t246 = -t251 * t272 + t255 * t256;
t262 = t246 * t260 + t259 * t264;
t261 = t259 * t270 + (t251 * t260 + t259 * t268) * t254;
t249 = qJD(4) + qJD(5);
t248 = cos(t250);
t242 = t261 * qJD(3);
t241 = t262 * qJD(3);
t240 = t263 * qJD(3);
t1 = [0, 0, 0, t241, t241 (t262 * t248 + (-t245 * t253 + t252 * t269) * t247) * t249 + (-t246 * t259 + t260 * t264) * t266; 0, 0, 0, t240, t240 (t263 * t248 + (-t243 * t253 - t256 * t269) * t247) * t249 + (-t244 * t259 + t260 * t265) * t266; 0, 0, 0, t242, t242 (t261 * t248 + (-t255 * t271 + t257 * t258) * t247) * t249 + (t260 * t270 + (-t251 * t259 + t260 * t268) * t254) * t266;];
JgD_rot  = t1;
