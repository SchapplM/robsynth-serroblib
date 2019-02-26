% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP7_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobia_rot_3_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:25
% EndTime: 2019-02-26 21:49:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (82->16), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->26)
t29 = cos(qJ(1));
t30 = t29 ^ 2;
t28 = cos(qJ(2));
t26 = sin(qJ(2));
t27 = sin(qJ(1));
t31 = t27 * t26;
t18 = atan2(-t31, -t28);
t16 = sin(t18);
t17 = cos(t18);
t14 = -t16 * t31 - t17 * t28;
t13 = 0.1e1 / t14 ^ 2;
t36 = t13 * t26;
t35 = t16 * t28;
t21 = t26 ^ 2;
t24 = 0.1e1 / t28 ^ 2;
t34 = t21 * t24;
t22 = t27 ^ 2;
t33 = t22 / t30;
t19 = 0.1e1 / (t22 * t34 + 0.1e1);
t32 = t27 * t19;
t23 = 0.1e1 / t28;
t20 = 0.1e1 / (t24 * t33 + 0.1e1);
t15 = (0.1e1 + t34) * t32;
t12 = 0.1e1 / t14;
t11 = 0.1e1 / (t30 * t21 * t13 + 0.1e1);
t1 = [t29 * t26 * t23 * t19, t15, 0, 0, 0, 0; (-t12 * t31 - (-t17 * t21 * t23 * t32 + (t19 - 0.1e1) * t26 * t16) * t30 * t36) * t11 (t28 * t12 - (-t27 * t35 + t17 * t26 + (-t17 * t31 + t35) * t15) * t36) * t29 * t11, 0, 0, 0, 0; (-0.1e1 - t33) * t23 * t20, -0.1e1 / t29 * t24 * t20 * t31, 0, 0, 0, 0;];
Ja_rot  = t1;
