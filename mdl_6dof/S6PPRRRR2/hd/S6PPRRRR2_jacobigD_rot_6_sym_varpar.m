% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR2
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
% Datum: 2019-02-26 19:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRR2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:16
% EndTime: 2019-02-26 19:43:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (72->23), mult. (242->56), div. (0->0), fcn. (284->12), ass. (0->32)
t241 = sin(pkin(12));
t247 = cos(pkin(6));
t263 = t241 * t247;
t242 = sin(pkin(7));
t243 = sin(pkin(6));
t262 = t242 * t243;
t261 = t242 * t247;
t246 = cos(pkin(7));
t260 = t243 * t246;
t244 = cos(pkin(13));
t259 = t244 * t246;
t245 = cos(pkin(12));
t258 = t245 * t247;
t248 = sin(qJ(4));
t257 = qJD(3) * t248;
t240 = sin(pkin(13));
t236 = -t241 * t240 + t244 * t258;
t256 = t236 * t246 - t245 * t262;
t238 = -t245 * t240 - t244 * t263;
t255 = t238 * t246 + t241 * t262;
t237 = t240 * t258 + t241 * t244;
t249 = sin(qJ(3));
t251 = cos(qJ(3));
t254 = t237 * t251 + t256 * t249;
t239 = -t240 * t263 + t245 * t244;
t253 = t239 * t251 + t255 * t249;
t252 = t249 * t261 + (t240 * t251 + t249 * t259) * t243;
t250 = cos(qJ(4));
t235 = (t252 * t250 + (-t244 * t262 + t247 * t246) * t248) * qJD(4) + (t251 * t261 + (-t240 * t249 + t251 * t259) * t243) * t257;
t234 = (t253 * t250 + (-t238 * t242 + t241 * t260) * t248) * qJD(4) + (-t239 * t249 + t255 * t251) * t257;
t233 = (t254 * t250 + (-t236 * t242 - t245 * t260) * t248) * qJD(4) + (-t237 * t249 + t256 * t251) * t257;
t1 = [0, 0, 0, t253 * qJD(3), t234, t234; 0, 0, 0, t254 * qJD(3), t233, t233; 0, 0, 0, t252 * qJD(3), t235, t235;];
JgD_rot  = t1;
