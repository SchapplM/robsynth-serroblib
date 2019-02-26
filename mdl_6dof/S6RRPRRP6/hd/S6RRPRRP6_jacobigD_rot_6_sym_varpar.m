% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRP6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:51
% EndTime: 2019-02-26 21:48:51
% DurationCPUTime: 0.09s
% Computational Cost: add. (45->22), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->24)
t231 = sin(pkin(11));
t233 = cos(pkin(11));
t236 = sin(qJ(2));
t239 = cos(qJ(2));
t242 = t236 * t231 - t239 * t233;
t249 = t242 * qJD(2);
t243 = t239 * t231 + t236 * t233;
t228 = t243 * qJD(2);
t232 = sin(pkin(6));
t237 = sin(qJ(1));
t248 = t232 * t237;
t240 = cos(qJ(1));
t247 = t232 * t240;
t246 = qJD(1) * t232;
t234 = cos(pkin(6));
t226 = t243 * t234;
t245 = t240 * t226 - t237 * t242;
t244 = -t237 * t226 - t240 * t242;
t238 = cos(qJ(4));
t235 = sin(qJ(4));
t225 = t242 * t234;
t224 = t234 * t249;
t223 = t234 * t228;
t1 = [0, t240 * t246, 0, -t237 * t223 - t240 * t249 + (-t225 * t240 - t237 * t243) * qJD(1) (t237 * t224 - t240 * t228) * t235 + (t235 * t248 + t244 * t238) * qJD(4) + (-t245 * t235 - t238 * t247) * qJD(1), 0; 0, t237 * t246, 0, t240 * t223 - t237 * t249 + (-t225 * t237 + t240 * t243) * qJD(1) (-t240 * t224 - t237 * t228) * t235 + (-t235 * t247 + t245 * t238) * qJD(4) + (t244 * t235 - t238 * t248) * qJD(1), 0; 0, 0, 0, t232 * t228, t234 * qJD(4) * t235 + (t243 * qJD(4) * t238 - t235 * t249) * t232, 0;];
JgD_rot  = t1;
