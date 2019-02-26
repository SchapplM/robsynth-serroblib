% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:47
% EndTime: 2019-02-26 21:55:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (77->23), mult. (197->50), div. (0->0), fcn. (216->10), ass. (0->29)
t242 = sin(pkin(12));
t244 = cos(pkin(12));
t246 = sin(qJ(2));
t248 = cos(qJ(2));
t251 = t246 * t242 - t248 * t244;
t259 = t251 * qJD(2);
t252 = t248 * t242 + t246 * t244;
t235 = t252 * qJD(2);
t241 = qJ(4) + qJ(5);
t239 = cos(t241);
t240 = qJD(4) + qJD(5);
t258 = t239 * t240;
t243 = sin(pkin(6));
t257 = t243 * t239;
t256 = qJD(1) * t243;
t245 = cos(pkin(6));
t255 = t240 * t243 + t245 * t259;
t233 = t252 * t245;
t247 = sin(qJ(1));
t249 = cos(qJ(1));
t254 = t249 * t233 - t247 * t251;
t253 = -t247 * t233 - t249 * t251;
t238 = sin(t241);
t232 = t251 * t245;
t230 = t245 * t235;
t229 = t243 * t235;
t228 = t249 * t230 - t247 * t259 + (-t232 * t247 + t249 * t252) * qJD(1);
t227 = -t247 * t230 - t249 * t259 + (-t232 * t249 - t247 * t252) * qJD(1);
t1 = [0, t249 * t256, 0, t227, t227, t253 * t258 + (-t249 * t235 + t255 * t247) * t238 + (-t254 * t238 - t249 * t257) * qJD(1); 0, t247 * t256, 0, t228, t228, t254 * t258 + (-t247 * t235 - t255 * t249) * t238 + (t253 * t238 - t247 * t257) * qJD(1); 0, 0, 0, t229, t229, t245 * t240 * t238 + (-t238 * t259 + t252 * t258) * t243;];
JgD_rot  = t1;
