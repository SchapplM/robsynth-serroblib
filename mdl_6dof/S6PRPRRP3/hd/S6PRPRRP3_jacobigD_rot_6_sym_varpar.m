% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JgD_rot = S6PRPRRP3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:34
% EndTime: 2019-02-26 19:51:34
% DurationCPUTime: 0.04s
% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
t143 = pkin(11) + qJ(4);
t141 = sin(t143);
t145 = sin(pkin(6));
t156 = t145 * t141;
t147 = cos(pkin(6));
t148 = sin(qJ(2));
t155 = t147 * t148;
t149 = cos(qJ(2));
t154 = t147 * t149;
t153 = qJD(2) * t141;
t152 = qJD(2) * t145;
t144 = sin(pkin(10));
t146 = cos(pkin(10));
t151 = t144 * t149 + t146 * t155;
t150 = -t144 * t155 + t146 * t149;
t142 = cos(t143);
t1 = [0, 0, 0, t150 * qJD(2) (t150 * t142 + t144 * t156) * qJD(4) + (-t144 * t154 - t146 * t148) * t153, 0; 0, 0, 0, t151 * qJD(2) (t151 * t142 - t146 * t156) * qJD(4) + (-t144 * t148 + t146 * t154) * t153, 0; 0, 0, 0, t148 * t152, t149 * t141 * t152 + (t142 * t145 * t148 + t141 * t147) * qJD(4), 0;];
JgD_rot  = t1;
