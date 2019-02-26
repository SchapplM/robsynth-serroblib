% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPPRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:35
% EndTime: 2019-02-26 20:22:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (172->19), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
t37 = qJ(1) + pkin(9);
t36 = cos(t37);
t55 = t36 ^ 2;
t42 = sin(qJ(5));
t35 = sin(t37);
t44 = cos(qJ(5));
t50 = t35 * t44;
t34 = atan2(t50, t42);
t31 = sin(t34);
t32 = cos(t34);
t26 = t31 * t50 + t32 * t42;
t25 = 0.1e1 / t26 ^ 2;
t54 = t25 * t44;
t41 = sin(qJ(6));
t43 = cos(qJ(6));
t46 = t42 * t43;
t30 = -t35 * t41 + t36 * t46;
t28 = 0.1e1 / t30 ^ 2;
t47 = t41 * t42;
t29 = t35 * t43 + t36 * t47;
t53 = t28 * t29;
t52 = t31 * t42;
t40 = t44 ^ 2;
t48 = 0.1e1 / t42 ^ 2 * t40;
t33 = 0.1e1 / (t35 ^ 2 * t48 + 0.1e1);
t51 = t35 * t33;
t49 = t36 * t44;
t45 = t29 ^ 2 * t28 + 0.1e1;
t38 = 0.1e1 / t42;
t27 = 0.1e1 / t30;
t24 = 0.1e1 / t26;
t23 = (-0.1e1 - t48) * t51;
t22 = 0.1e1 / t45;
t21 = 0.1e1 / (t55 * t40 * t25 + 0.1e1);
t1 = [t38 * t33 * t49, 0, 0, 0, t23, 0; (t24 * t50 + (t32 * t38 * t40 * t51 + (-t33 + 0.1e1) * t44 * t31) * t55 * t54) * t21, 0, 0, 0 (t42 * t24 + (-t35 * t52 + t32 * t44 + (t32 * t50 - t52) * t23) * t54) * t36 * t21, 0; ((-t35 * t47 + t36 * t43) * t27 - (-t35 * t46 - t36 * t41) * t53) * t22, 0, 0, 0 (t27 * t41 - t43 * t53) * t22 * t49, t45 * t22;];
Ja_rot  = t1;
