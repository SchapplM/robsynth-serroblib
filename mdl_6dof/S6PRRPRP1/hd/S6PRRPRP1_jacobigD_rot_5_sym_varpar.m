% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRP1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobigD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:21
% EndTime: 2019-02-26 20:01:21
% DurationCPUTime: 0.04s
% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
t139 = qJ(3) + pkin(11);
t137 = sin(t139);
t141 = sin(pkin(6));
t152 = t141 * t137;
t143 = cos(pkin(6));
t144 = sin(qJ(2));
t151 = t143 * t144;
t145 = cos(qJ(2));
t150 = t143 * t145;
t149 = qJD(2) * t137;
t148 = qJD(2) * t141;
t140 = sin(pkin(10));
t142 = cos(pkin(10));
t147 = t140 * t145 + t142 * t151;
t146 = -t140 * t151 + t142 * t145;
t138 = cos(t139);
t1 = [0, 0, t146 * qJD(2), 0 (t146 * t138 + t140 * t152) * qJD(3) + (-t140 * t150 - t142 * t144) * t149, 0; 0, 0, t147 * qJD(2), 0 (t147 * t138 - t142 * t152) * qJD(3) + (-t140 * t144 + t142 * t150) * t149, 0; 0, 0, t144 * t148, 0, t145 * t137 * t148 + (t138 * t141 * t144 + t137 * t143) * qJD(3), 0;];
JgD_rot  = t1;
