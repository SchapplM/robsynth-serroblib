% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:16
% EndTime: 2019-02-26 20:26:16
% DurationCPUTime: 0.06s
% Computational Cost: add. (263->18), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->29)
t35 = qJ(1) + pkin(9);
t33 = cos(t35);
t29 = t33 ^ 2;
t34 = pkin(10) + qJ(4);
t32 = cos(t34);
t30 = sin(t34);
t31 = sin(t35);
t38 = t31 * t30;
t22 = atan2(-t38, -t32);
t20 = sin(t22);
t21 = cos(t22);
t18 = -t20 * t38 - t21 * t32;
t17 = 0.1e1 / t18 ^ 2;
t44 = t17 * t30;
t43 = t20 * t32;
t25 = t30 ^ 2;
t37 = t32 ^ 2;
t42 = t25 / t37;
t36 = t31 ^ 2;
t41 = 0.1e1 / t36 * t29;
t40 = t30 * t33;
t23 = 0.1e1 / (t36 * t42 + 0.1e1);
t39 = t31 * t23;
t27 = 0.1e1 / t32;
t24 = 0.1e1 / (t37 * t41 + 0.1e1);
t19 = (0.1e1 + t42) * t39;
t16 = 0.1e1 / t18;
t15 = 0.1e1 / (t29 * t25 * t17 + 0.1e1);
t1 = [t27 * t23 * t40, 0, 0, t19, 0, 0; (-t16 * t38 - (-t21 * t25 * t27 * t39 + (t23 - 0.1e1) * t30 * t20) * t29 * t44) * t15, 0, 0 (t32 * t16 - (-t31 * t43 + t21 * t30 + (-t21 * t38 + t43) * t19) * t44) * t33 * t15, 0, 0; (-0.1e1 - t41) * t32 * t24, 0, 0, -0.1e1 / t31 * t24 * t40, 0, 0;];
Ja_rot  = t1;
