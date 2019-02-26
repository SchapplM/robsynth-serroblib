% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR7_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:05
% EndTime: 2019-02-26 20:14:05
% DurationCPUTime: 0.08s
% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
t230 = sin(pkin(7));
t231 = sin(pkin(6));
t249 = t230 * t231;
t233 = cos(pkin(7));
t237 = cos(qJ(3));
t248 = t233 * t237;
t234 = cos(pkin(6));
t236 = sin(qJ(2));
t247 = t234 * t236;
t238 = cos(qJ(2));
t246 = t234 * t238;
t235 = sin(qJ(3));
t245 = t235 * t238;
t244 = t236 * t237;
t243 = qJD(2) * t235;
t229 = sin(pkin(12));
t232 = cos(pkin(12));
t242 = -t229 * t236 + t232 * t246;
t241 = t229 * t238 + t232 * t247;
t240 = -t229 * t246 - t232 * t236;
t239 = t229 * t247 - t232 * t238;
t228 = t239 * qJD(2);
t227 = t241 * qJD(2);
t1 = [0, 0, -t228 * t230, -t228 * t248 + t240 * t243 + (-t239 * t237 + (t229 * t249 + t240 * t233) * t235) * qJD(3), 0, 0; 0, 0, t227 * t230, t227 * t248 + t242 * t243 + (t241 * t237 + (-t232 * t249 + t242 * t233) * t235) * qJD(3), 0, 0; 0, 0, qJD(2) * t236 * t249, t234 * t230 * qJD(3) * t235 + ((t233 * t245 + t244) * qJD(3) + (t233 * t244 + t245) * qJD(2)) * t231, 0, 0;];
JgD_rot  = t1;
