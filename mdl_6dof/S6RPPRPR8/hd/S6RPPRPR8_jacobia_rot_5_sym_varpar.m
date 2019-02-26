% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR8_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:44
% EndTime: 2019-02-26 20:29:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (168->15), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->28)
t32 = sin(qJ(1));
t30 = t32 ^ 2;
t29 = pkin(9) + qJ(4);
t27 = sin(t29);
t28 = cos(t29);
t33 = cos(qJ(1));
t36 = t33 * t28;
t21 = atan2(-t36, t27);
t19 = sin(t21);
t20 = cos(t21);
t17 = -t19 * t36 + t20 * t27;
t16 = 0.1e1 / t17 ^ 2;
t42 = t16 * t28;
t41 = t19 * t27;
t26 = t28 ^ 2;
t34 = t27 ^ 2;
t40 = 0.1e1 / t34 * t26;
t39 = t28 * t32;
t35 = t33 ^ 2;
t38 = t30 / t35;
t22 = 0.1e1 / (t35 * t40 + 0.1e1);
t37 = t33 * t22;
t24 = 0.1e1 / t27;
t23 = 0.1e1 / (t34 * t38 + 0.1e1);
t18 = (0.1e1 + t40) * t37;
t15 = 0.1e1 / t17;
t14 = 0.1e1 / (t30 * t26 * t16 + 0.1e1);
t1 = [t24 * t22 * t39, 0, 0, t18, 0, 0; (-t15 * t36 + (-t20 * t24 * t26 * t37 + (-t22 + 0.1e1) * t28 * t19) * t30 * t42) * t14, 0, 0 (t27 * t15 + (t33 * t41 + t20 * t28 + (-t20 * t36 - t41) * t18) * t42) * t32 * t14, 0, 0; (0.1e1 + t38) * t27 * t23, 0, 0, 0.1e1 / t33 * t23 * t39, 0, 0;];
Ja_rot  = t1;
