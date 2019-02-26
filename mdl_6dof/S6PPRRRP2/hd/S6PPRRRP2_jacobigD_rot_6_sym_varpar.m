% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRP2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:05
% EndTime: 2019-02-26 19:42:05
% DurationCPUTime: 0.11s
% Computational Cost: add. (41->23), mult. (141->56), div. (0->0), fcn. (164->12), ass. (0->29)
t236 = sin(pkin(11));
t242 = cos(pkin(6));
t258 = t236 * t242;
t237 = sin(pkin(7));
t238 = sin(pkin(6));
t257 = t237 * t238;
t256 = t237 * t242;
t241 = cos(pkin(7));
t255 = t238 * t241;
t239 = cos(pkin(12));
t254 = t239 * t241;
t240 = cos(pkin(11));
t253 = t240 * t242;
t243 = sin(qJ(4));
t252 = qJD(3) * t243;
t235 = sin(pkin(12));
t231 = -t235 * t236 + t239 * t253;
t251 = t231 * t241 - t240 * t257;
t233 = -t235 * t240 - t239 * t258;
t250 = t233 * t241 + t236 * t257;
t232 = t235 * t253 + t236 * t239;
t244 = sin(qJ(3));
t246 = cos(qJ(3));
t249 = t232 * t246 + t251 * t244;
t234 = -t235 * t258 + t239 * t240;
t248 = t234 * t246 + t250 * t244;
t247 = t244 * t256 + (t235 * t246 + t244 * t254) * t238;
t245 = cos(qJ(4));
t1 = [0, 0, 0, t248 * qJD(3) (t248 * t245 + (-t233 * t237 + t236 * t255) * t243) * qJD(4) + (-t234 * t244 + t250 * t246) * t252, 0; 0, 0, 0, t249 * qJD(3) (t249 * t245 + (-t231 * t237 - t240 * t255) * t243) * qJD(4) + (-t232 * t244 + t251 * t246) * t252, 0; 0, 0, 0, t247 * qJD(3) (t247 * t245 + (-t239 * t257 + t241 * t242) * t243) * qJD(4) + (t246 * t256 + (-t235 * t244 + t246 * t254) * t238) * t252, 0;];
JgD_rot  = t1;
