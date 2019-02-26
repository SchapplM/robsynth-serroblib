% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPPR2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:46
% EndTime: 2019-02-26 19:58:46
% DurationCPUTime: 0.04s
% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
t141 = qJ(3) + pkin(11);
t140 = cos(t141);
t143 = sin(pkin(6));
t154 = t143 * t140;
t145 = cos(pkin(6));
t146 = sin(qJ(2));
t153 = t145 * t146;
t147 = cos(qJ(2));
t152 = t145 * t147;
t151 = qJD(2) * t140;
t150 = qJD(2) * t143;
t142 = sin(pkin(10));
t144 = cos(pkin(10));
t149 = t142 * t147 + t144 * t153;
t148 = -t142 * t153 + t144 * t147;
t139 = sin(t141);
t1 = [0, 0, t148 * qJD(2), 0, 0 (-t148 * t139 + t142 * t154) * qJD(3) + (-t142 * t152 - t144 * t146) * t151; 0, 0, t149 * qJD(2), 0, 0 (-t149 * t139 - t144 * t154) * qJD(3) + (-t142 * t146 + t144 * t152) * t151; 0, 0, t146 * t150, 0, 0, t147 * t140 * t150 + (-t139 * t143 * t146 + t140 * t145) * qJD(3);];
JgD_rot  = t1;
