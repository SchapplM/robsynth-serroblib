% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14V3_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobia_rot_3_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:04
% DurationCPUTime: 0.06s
% Computational Cost: add. (82->16), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->26)
t27 = cos(qJ(1));
t28 = t27 ^ 2;
t26 = cos(qJ(2));
t24 = sin(qJ(2));
t25 = sin(qJ(1));
t29 = t25 * t24;
t16 = atan2(-t29, -t26);
t14 = sin(t16);
t15 = cos(t16);
t12 = -t14 * t29 - t15 * t26;
t11 = 0.1e1 / t12 ^ 2;
t34 = t11 * t24;
t33 = t14 * t26;
t19 = t24 ^ 2;
t22 = 0.1e1 / t26 ^ 2;
t32 = t19 * t22;
t20 = t25 ^ 2;
t31 = t20 / t28;
t17 = 0.1e1 / (t20 * t32 + 0.1e1);
t30 = t25 * t17;
t21 = 0.1e1 / t26;
t18 = 0.1e1 / (t22 * t31 + 0.1e1);
t13 = (0.1e1 + t32) * t30;
t10 = 0.1e1 / t12;
t9 = 0.1e1 / (t28 * t19 * t11 + 0.1e1);
t1 = [t27 * t24 * t21 * t17, t13, 0, 0, 0, 0; (-t10 * t29 - (-t15 * t19 * t21 * t30 + (t17 - 0.1e1) * t24 * t14) * t28 * t34) * t9 (t26 * t10 - (-t25 * t33 + t15 * t24 + (-t15 * t29 + t33) * t13) * t34) * t9 * t27, 0, 0, 0, 0; (-0.1e1 - t31) * t21 * t18, -0.1e1 / t27 * t22 * t18 * t29, 0, 0, 0, 0;];
Ja_rot  = t1;
