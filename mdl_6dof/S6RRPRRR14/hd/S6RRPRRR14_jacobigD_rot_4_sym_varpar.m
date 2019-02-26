% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR14_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobigD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:45
% EndTime: 2019-02-26 22:55:46
% DurationCPUTime: 0.07s
% Computational Cost: add. (28->18), mult. (102->43), div. (0->0), fcn. (102->12), ass. (0->29)
t241 = cos(pkin(14));
t243 = cos(pkin(7));
t262 = t243 * t241;
t245 = sin(qJ(2));
t246 = sin(qJ(1));
t261 = t245 * t246;
t248 = cos(qJ(1));
t260 = t245 * t248;
t247 = cos(qJ(2));
t259 = t246 * t247;
t258 = t247 * t248;
t240 = sin(pkin(6));
t257 = qJD(1) * t240;
t237 = sin(pkin(14));
t256 = qJD(2) * t237;
t239 = sin(pkin(7));
t255 = t240 * t239 * t241;
t254 = t246 * t257;
t253 = t248 * t257;
t244 = cos(pkin(6));
t252 = t244 * t258 - t261;
t251 = -t244 * t259 - t260;
t250 = -t244 * t260 - t259;
t249 = t244 * t261 - t258;
t242 = cos(pkin(8));
t238 = sin(pkin(8));
t236 = t251 * qJD(1) + t250 * qJD(2);
t235 = -t252 * qJD(1) + t249 * qJD(2);
t1 = [0, t253, 0 -(t235 * t262 - t251 * t256 + (-t250 * t237 + t248 * t255) * qJD(1)) * t238 + (-t235 * t239 + t243 * t253) * t242, 0, 0; 0, t254, 0 -(t236 * t262 - t252 * t256 + (t249 * t237 + t246 * t255) * qJD(1)) * t238 + (-t236 * t239 + t243 * t254) * t242, 0, 0; 0, 0, 0 (-(-t237 * t247 - t245 * t262) * t238 + t245 * t239 * t242) * t240 * qJD(2), 0, 0;];
JgD_rot  = t1;
