% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRP4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:04
% EndTime: 2019-02-26 20:03:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t137 = sin(pkin(6));
t141 = sin(qJ(2));
t150 = t137 * t141;
t142 = cos(qJ(3));
t149 = t137 * t142;
t139 = cos(pkin(6));
t148 = t139 * t141;
t143 = cos(qJ(2));
t147 = t139 * t143;
t146 = qJD(2) * t142;
t136 = sin(pkin(10));
t138 = cos(pkin(10));
t145 = t136 * t143 + t138 * t148;
t144 = -t136 * t148 + t138 * t143;
t140 = sin(qJ(3));
t1 = [0, 0, t144 * qJD(2), 0 (t136 * t149 - t144 * t140) * qJD(3) + (-t136 * t147 - t138 * t141) * t146, 0; 0, 0, t145 * qJD(2), 0 (-t138 * t149 - t145 * t140) * qJD(3) + (-t136 * t141 + t138 * t147) * t146, 0; 0, 0, qJD(2) * t150, 0, t137 * t143 * t146 + (t139 * t142 - t140 * t150) * qJD(3), 0;];
JgD_rot  = t1;
