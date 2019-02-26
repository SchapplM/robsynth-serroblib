% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR1_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:06
% EndTime: 2019-02-26 21:28:06
% DurationCPUTime: 0.06s
% Computational Cost: add. (193->17), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->27)
t36 = cos(qJ(1));
t37 = t36 ^ 2;
t32 = qJ(2) + pkin(10);
t31 = cos(t32);
t30 = sin(t32);
t35 = sin(qJ(1));
t38 = t35 * t30;
t24 = atan2(-t38, -t31);
t22 = sin(t24);
t23 = cos(t24);
t20 = -t22 * t38 - t23 * t31;
t19 = 0.1e1 / t20 ^ 2;
t43 = t19 * t30;
t42 = t22 * t31;
t27 = t30 ^ 2;
t29 = 0.1e1 / t31 ^ 2;
t41 = t27 * t29;
t33 = t35 ^ 2;
t40 = t33 / t37;
t25 = 0.1e1 / (t33 * t41 + 0.1e1);
t39 = t35 * t25;
t28 = 0.1e1 / t31;
t26 = 0.1e1 / (t29 * t40 + 0.1e1);
t21 = (0.1e1 + t41) * t39;
t18 = 0.1e1 / t20;
t17 = 0.1e1 / (t37 * t27 * t19 + 0.1e1);
t1 = [t36 * t30 * t28 * t25, t21, 0, 0, 0, 0; (-t18 * t38 - (-t23 * t27 * t28 * t39 + (t25 - 0.1e1) * t30 * t22) * t37 * t43) * t17 (t31 * t18 - (-t35 * t42 + t23 * t30 + (-t23 * t38 + t42) * t21) * t43) * t36 * t17, 0, 0, 0, 0; (-0.1e1 - t40) * t28 * t26, -0.1e1 / t36 * t29 * t26 * t38, 0, 0, 0, 0;];
Ja_rot  = t1;
