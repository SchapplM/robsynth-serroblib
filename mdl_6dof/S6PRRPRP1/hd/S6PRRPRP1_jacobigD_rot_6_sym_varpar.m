% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JgD_rot = S6PRRPRP1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:11
% EndTime: 2019-02-26 20:01:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
t146 = qJ(3) + pkin(11);
t144 = sin(t146);
t148 = sin(pkin(6));
t159 = t148 * t144;
t150 = cos(pkin(6));
t151 = sin(qJ(2));
t158 = t150 * t151;
t152 = cos(qJ(2));
t157 = t150 * t152;
t156 = qJD(2) * t144;
t155 = qJD(2) * t148;
t147 = sin(pkin(10));
t149 = cos(pkin(10));
t154 = t147 * t152 + t149 * t158;
t153 = -t147 * t158 + t149 * t152;
t145 = cos(t146);
t1 = [0, 0, t153 * qJD(2), 0 (t153 * t145 + t147 * t159) * qJD(3) + (-t147 * t157 - t149 * t151) * t156, 0; 0, 0, t154 * qJD(2), 0 (t154 * t145 - t149 * t159) * qJD(3) + (-t147 * t151 + t149 * t157) * t156, 0; 0, 0, t151 * t155, 0, t152 * t144 * t155 + (t145 * t148 * t151 + t144 * t150) * qJD(3), 0;];
JgD_rot  = t1;
