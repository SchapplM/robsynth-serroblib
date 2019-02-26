% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRP8_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:49
% EndTime: 2019-02-26 22:43:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (45->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
t240 = sin(pkin(6));
t243 = sin(qJ(1));
t258 = t240 * t243;
t245 = cos(qJ(1));
t257 = t240 * t245;
t242 = sin(qJ(2));
t256 = t242 * t243;
t255 = t242 * t245;
t244 = cos(qJ(2));
t254 = t243 * t244;
t253 = t245 * t244;
t252 = qJD(1) * t240;
t239 = qJ(3) + qJ(4);
t236 = sin(t239);
t251 = qJD(2) * t236;
t250 = qJD(2) * t240;
t241 = cos(pkin(6));
t249 = t241 * t253 - t256;
t248 = t241 * t254 + t255;
t247 = t241 * t255 + t254;
t246 = -t241 * t256 + t253;
t238 = qJD(3) + qJD(4);
t237 = cos(t239);
t235 = t242 * t250;
t234 = t248 * qJD(1) + t247 * qJD(2);
t233 = t249 * qJD(1) + t246 * qJD(2);
t1 = [0, t245 * t252, t233, t233 (t236 * t258 + t246 * t237) * t238 - t248 * t251 + (-t247 * t236 - t237 * t257) * qJD(1), 0; 0, t243 * t252, t234, t234 (-t236 * t257 + t247 * t237) * t238 + t249 * t251 + (t246 * t236 - t237 * t258) * qJD(1), 0; 0, 0, t235, t235, t240 * t242 * t238 * t237 + (t238 * t241 + t244 * t250) * t236, 0;];
JgD_rot  = t1;
