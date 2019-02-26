% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPPR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:13
% EndTime: 2019-02-26 19:58:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
t145 = qJ(3) + pkin(11);
t143 = sin(t145);
t147 = sin(pkin(6));
t158 = t147 * t143;
t149 = cos(pkin(6));
t150 = sin(qJ(2));
t157 = t149 * t150;
t151 = cos(qJ(2));
t156 = t149 * t151;
t155 = qJD(2) * t143;
t154 = qJD(2) * t147;
t146 = sin(pkin(10));
t148 = cos(pkin(10));
t153 = t146 * t151 + t148 * t157;
t152 = -t146 * t157 + t148 * t151;
t144 = cos(t145);
t1 = [0, 0, t152 * qJD(2), 0, 0 (t152 * t144 + t146 * t158) * qJD(3) + (-t146 * t156 - t148 * t150) * t155; 0, 0, t153 * qJD(2), 0, 0 (t153 * t144 - t148 * t158) * qJD(3) + (-t146 * t150 + t148 * t156) * t155; 0, 0, t150 * t154, 0, 0, t151 * t143 * t154 + (t144 * t147 * t150 + t143 * t149) * qJD(3);];
JgD_rot  = t1;
