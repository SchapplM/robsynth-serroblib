% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRP3_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobigD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:34
% EndTime: 2019-02-26 19:51:34
% DurationCPUTime: 0.06s
% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
t135 = pkin(11) + qJ(4);
t133 = sin(t135);
t137 = sin(pkin(6));
t148 = t137 * t133;
t139 = cos(pkin(6));
t140 = sin(qJ(2));
t147 = t139 * t140;
t141 = cos(qJ(2));
t146 = t139 * t141;
t145 = qJD(2) * t133;
t144 = qJD(2) * t137;
t136 = sin(pkin(10));
t138 = cos(pkin(10));
t143 = t136 * t141 + t138 * t147;
t142 = -t136 * t147 + t138 * t141;
t134 = cos(t135);
t1 = [0, 0, 0, t142 * qJD(2) (t134 * t142 + t136 * t148) * qJD(4) + (-t136 * t146 - t138 * t140) * t145, 0; 0, 0, 0, t143 * qJD(2) (t134 * t143 - t138 * t148) * qJD(4) + (-t136 * t140 + t138 * t146) * t145, 0; 0, 0, 0, t140 * t144, t141 * t133 * t144 + (t134 * t137 * t140 + t133 * t139) * qJD(4), 0;];
JgD_rot  = t1;
