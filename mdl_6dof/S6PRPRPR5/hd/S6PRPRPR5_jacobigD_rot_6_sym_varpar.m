% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR5_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:45
% EndTime: 2019-02-26 19:48:46
% DurationCPUTime: 0.07s
% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
t137 = pkin(11) + qJ(4);
t136 = cos(t137);
t139 = sin(pkin(6));
t150 = t139 * t136;
t141 = cos(pkin(6));
t142 = sin(qJ(2));
t149 = t141 * t142;
t143 = cos(qJ(2));
t148 = t141 * t143;
t147 = qJD(2) * t136;
t146 = qJD(2) * t139;
t138 = sin(pkin(10));
t140 = cos(pkin(10));
t145 = t138 * t143 + t140 * t149;
t144 = -t138 * t149 + t140 * t143;
t135 = sin(t137);
t1 = [0, 0, 0, t144 * qJD(2), 0 (-t135 * t144 + t138 * t150) * qJD(4) + (-t138 * t148 - t140 * t142) * t147; 0, 0, 0, t145 * qJD(2), 0 (-t135 * t145 - t140 * t150) * qJD(4) + (-t138 * t142 + t140 * t148) * t147; 0, 0, 0, t142 * t146, 0, t143 * t136 * t146 + (-t135 * t139 * t142 + t136 * t141) * qJD(4);];
JgD_rot  = t1;
