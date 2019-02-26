% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:15
% EndTime: 2019-02-26 19:48:15
% DurationCPUTime: 0.04s
% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
t142 = pkin(11) + qJ(4);
t140 = sin(t142);
t144 = sin(pkin(6));
t155 = t144 * t140;
t146 = cos(pkin(6));
t147 = sin(qJ(2));
t154 = t146 * t147;
t148 = cos(qJ(2));
t153 = t146 * t148;
t152 = qJD(2) * t140;
t151 = qJD(2) * t144;
t143 = sin(pkin(10));
t145 = cos(pkin(10));
t150 = t143 * t148 + t145 * t154;
t149 = -t143 * t154 + t145 * t148;
t141 = cos(t142);
t1 = [0, 0, 0, t149 * qJD(2), 0 (t149 * t141 + t143 * t155) * qJD(4) + (-t143 * t153 - t145 * t147) * t152; 0, 0, 0, t150 * qJD(2), 0 (t150 * t141 - t145 * t155) * qJD(4) + (-t143 * t147 + t145 * t153) * t152; 0, 0, 0, t147 * t151, 0, t148 * t140 * t151 + (t141 * t144 * t147 + t140 * t146) * qJD(4);];
JgD_rot  = t1;
