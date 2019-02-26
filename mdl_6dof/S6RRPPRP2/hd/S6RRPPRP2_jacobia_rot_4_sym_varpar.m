% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP2_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobia_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:21
% EndTime: 2019-02-26 21:25:21
% DurationCPUTime: 0.06s
% Computational Cost: add. (192->17), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->28)
t36 = cos(qJ(1));
t34 = t36 ^ 2;
t32 = qJ(2) + pkin(9);
t31 = cos(t32);
t30 = sin(t32);
t35 = sin(qJ(1));
t39 = t35 * t30;
t24 = atan2(-t39, -t31);
t22 = sin(t24);
t23 = cos(t24);
t20 = -t22 * t39 - t23 * t31;
t19 = 0.1e1 / t20 ^ 2;
t45 = t19 * t30;
t44 = t22 * t31;
t27 = t30 ^ 2;
t37 = t31 ^ 2;
t43 = t27 / t37;
t42 = t30 * t36;
t38 = t35 ^ 2;
t41 = 0.1e1 / t38 * t34;
t25 = 0.1e1 / (t38 * t43 + 0.1e1);
t40 = t35 * t25;
t28 = 0.1e1 / t31;
t26 = 0.1e1 / (t37 * t41 + 0.1e1);
t21 = (0.1e1 + t43) * t40;
t18 = 0.1e1 / t20;
t17 = 0.1e1 / (t34 * t27 * t19 + 0.1e1);
t1 = [t28 * t25 * t42, t21, 0, 0, 0, 0; (-t18 * t39 - (-t23 * t27 * t28 * t40 + (t25 - 0.1e1) * t30 * t22) * t34 * t45) * t17 (t31 * t18 - (-t35 * t44 + t23 * t30 + (-t23 * t39 + t44) * t21) * t45) * t36 * t17, 0, 0, 0, 0; (-0.1e1 - t41) * t31 * t26, -0.1e1 / t35 * t26 * t42, 0, 0, 0, 0;];
Ja_rot  = t1;
