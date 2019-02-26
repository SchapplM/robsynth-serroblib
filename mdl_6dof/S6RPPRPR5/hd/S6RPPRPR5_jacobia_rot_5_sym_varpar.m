% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:00
% EndTime: 2019-02-26 20:28:01
% DurationCPUTime: 0.07s
% Computational Cost: add. (53->18), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->34)
t33 = sin(qJ(4));
t34 = sin(qJ(1));
t35 = cos(qJ(4));
t40 = t34 * t35;
t27 = atan2(t40, t33);
t24 = sin(t27);
t25 = cos(t27);
t18 = t24 * t40 + t25 * t33;
t17 = 0.1e1 / t18 ^ 2;
t36 = cos(qJ(1));
t47 = t17 * t36 ^ 2;
t32 = cos(pkin(9));
t37 = t36 * t32;
t31 = sin(pkin(9));
t42 = t34 * t31;
t23 = t33 * t37 - t42;
t21 = 0.1e1 / t23 ^ 2;
t38 = t36 * t31;
t41 = t34 * t32;
t22 = t33 * t38 + t41;
t46 = t21 * t22;
t45 = t24 * t33;
t30 = t35 ^ 2;
t44 = 0.1e1 / t33 ^ 2 * t30;
t26 = 0.1e1 / (t34 ^ 2 * t44 + 0.1e1);
t43 = t34 * t26;
t39 = t35 * t36;
t28 = 0.1e1 / t33;
t20 = 0.1e1 / t23;
t19 = (-0.1e1 - t44) * t43;
t16 = 0.1e1 / t18;
t15 = 0.1e1 / (t22 ^ 2 * t21 + 0.1e1);
t14 = 0.1e1 / (t30 * t47 + 0.1e1);
t1 = [t28 * t26 * t39, 0, 0, t19, 0, 0; (t16 * t40 + (t25 * t28 * t30 * t43 + (-t26 + 0.1e1) * t35 * t24) * t35 * t47) * t14, 0, 0 (t33 * t16 + (-t34 * t45 + t25 * t35 + (t25 * t40 - t45) * t19) * t35 * t17) * t36 * t14, 0, 0; ((-t33 * t42 + t37) * t20 - (-t33 * t41 - t38) * t46) * t15, 0, 0 (t20 * t31 - t32 * t46) * t15 * t39, 0, 0;];
Ja_rot  = t1;
