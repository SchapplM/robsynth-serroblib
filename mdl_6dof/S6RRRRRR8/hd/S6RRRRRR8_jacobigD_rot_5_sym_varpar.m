% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR8_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:33
% EndTime: 2019-02-26 22:51:33
% DurationCPUTime: 0.14s
% Computational Cost: add. (68->24), mult. (237->54), div. (0->0), fcn. (245->10), ass. (0->34)
t239 = sin(pkin(7));
t240 = sin(pkin(6));
t266 = t240 * t239;
t241 = cos(pkin(7));
t246 = cos(qJ(3));
t265 = t241 * t246;
t243 = sin(qJ(3));
t247 = cos(qJ(2));
t264 = t243 * t247;
t244 = sin(qJ(2));
t263 = t244 * t246;
t245 = sin(qJ(1));
t262 = t245 * t244;
t261 = t245 * t247;
t248 = cos(qJ(1));
t260 = t248 * t244;
t259 = t248 * t247;
t258 = qJD(1) * t240;
t257 = qJD(2) * t243;
t256 = t245 * t266;
t255 = t248 * t266;
t254 = t245 * t258;
t253 = t248 * t258;
t242 = cos(pkin(6));
t252 = t242 * t259 - t262;
t251 = -t242 * t261 - t260;
t250 = t242 * t260 + t261;
t249 = t242 * t262 - t259;
t238 = t251 * qJD(1) - t250 * qJD(2);
t237 = -t252 * qJD(1) + t249 * qJD(2);
t236 = t242 * t239 * qJD(3) * t243 + ((t241 * t264 + t263) * qJD(3) + (t241 * t263 + t264) * qJD(2)) * t240;
t235 = -t238 * t265 + t252 * t257 + (-t249 * t243 - t246 * t256) * qJD(1) + (t250 * t246 + (t252 * t241 - t255) * t243) * qJD(3);
t234 = -t237 * t265 + t251 * t257 + (-t250 * t243 - t246 * t255) * qJD(1) + (-t249 * t246 + (t251 * t241 + t256) * t243) * qJD(3);
t1 = [0, t253, -t237 * t239 + t241 * t253, t234, t234, 0; 0, t254, -t238 * t239 + t241 * t254, t235, t235, 0; 0, 0, qJD(2) * t244 * t266, t236, t236, 0;];
JgD_rot  = t1;
