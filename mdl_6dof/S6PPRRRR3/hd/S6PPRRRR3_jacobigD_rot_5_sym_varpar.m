% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRR3_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobigD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:56
% EndTime: 2019-02-26 19:43:56
% DurationCPUTime: 0.15s
% Computational Cost: add. (70->30), mult. (233->69), div. (0->0), fcn. (270->14), ass. (0->38)
t242 = sin(pkin(14));
t246 = sin(pkin(6));
t253 = sin(qJ(3));
t255 = cos(qJ(3));
t247 = cos(pkin(14));
t250 = cos(pkin(7));
t266 = t247 * t250;
t245 = sin(pkin(7));
t251 = cos(pkin(6));
t268 = t245 * t251;
t275 = (t242 * t255 + t253 * t266) * t246 + t253 * t268;
t248 = cos(pkin(13));
t243 = sin(pkin(13));
t270 = t243 * t251;
t241 = -t242 * t270 + t248 * t247;
t240 = -t248 * t242 - t247 * t270;
t269 = t245 * t246;
t260 = t240 * t250 + t243 * t269;
t274 = t241 * t255 + t260 * t253;
t265 = t248 * t251;
t239 = t242 * t265 + t243 * t247;
t238 = -t243 * t242 + t247 * t265;
t261 = -t238 * t250 + t248 * t269;
t273 = -t239 * t255 + t261 * t253;
t267 = t246 * t250;
t249 = cos(pkin(8));
t254 = cos(qJ(4));
t264 = t249 * t254;
t252 = sin(qJ(4));
t263 = qJD(3) * t252;
t258 = -t239 * t253 - t261 * t255;
t257 = -t241 * t253 + t260 * t255;
t256 = t255 * t268 + (-t242 * t253 + t255 * t266) * t246;
t244 = sin(pkin(8));
t237 = t275 * qJD(3);
t236 = t274 * qJD(3);
t235 = t273 * qJD(3);
t1 = [0, 0, 0, t236 * t244, t236 * t264 + (t274 * t254 + (t257 * t249 + (-t240 * t245 + t243 * t267) * t244) * t252) * qJD(4) + t257 * t263, 0; 0, 0, 0, -t235 * t244, -t235 * t264 + (-t273 * t254 + (t258 * t249 + (-t238 * t245 - t248 * t267) * t244) * t252) * qJD(4) + t258 * t263, 0; 0, 0, 0, t237 * t244, t237 * t264 + (t275 * t254 + (t256 * t249 + (-t247 * t269 + t251 * t250) * t244) * t252) * qJD(4) + t256 * t263, 0;];
JgD_rot  = t1;
