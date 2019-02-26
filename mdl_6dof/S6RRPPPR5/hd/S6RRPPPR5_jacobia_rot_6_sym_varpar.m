% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPPR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:18
% EndTime: 2019-02-26 21:24:18
% DurationCPUTime: 0.11s
% Computational Cost: add. (157->26), mult. (347->67), div. (52->9), fcn. (503->11), ass. (0->39)
t50 = cos(pkin(9));
t53 = sin(qJ(1));
t55 = cos(qJ(2));
t49 = sin(pkin(9));
t56 = cos(qJ(1));
t59 = t56 * t49;
t40 = -t53 * t50 + t55 * t59;
t58 = t56 * t50;
t41 = t53 * t49 + t55 * t58;
t51 = sin(qJ(6));
t54 = cos(qJ(6));
t33 = t40 * t54 + t41 * t51;
t31 = 0.1e1 / t33 ^ 2;
t32 = t40 * t51 - t41 * t54;
t66 = t31 * t32;
t52 = sin(qJ(2));
t61 = t53 * t52;
t44 = atan2(-t61, -t55);
t42 = sin(t44);
t43 = cos(t44);
t36 = -t42 * t61 - t43 * t55;
t35 = 0.1e1 / t36 ^ 2;
t65 = t35 * t56 ^ 2;
t46 = t52 ^ 2;
t64 = t46 / t55 ^ 2;
t63 = t52 * t56;
t45 = 0.1e1 / (t53 ^ 2 * t64 + 0.1e1);
t62 = t53 * t45;
t60 = t53 * t55;
t57 = t32 ^ 2 * t31 + 0.1e1;
t47 = 0.1e1 / t55;
t39 = -t50 * t60 + t59;
t38 = -t49 * t60 - t58;
t37 = (0.1e1 + t64) * t62;
t34 = 0.1e1 / t36;
t30 = 0.1e1 / t33;
t29 = 0.1e1 / (t46 * t65 + 0.1e1);
t28 = 0.1e1 / t57;
t1 = [t47 * t45 * t63, t37, 0, 0, 0, 0; (-t34 * t61 - (-t43 * t46 * t47 * t62 + (t45 - 0.1e1) * t52 * t42) * t52 * t65) * t29 (t55 * t34 - (-t42 * t60 + t43 * t52 + (t42 * t55 - t43 * t61) * t37) * t52 * t35) * t56 * t29, 0, 0, 0, 0; ((t38 * t51 - t39 * t54) * t30 - (t38 * t54 + t39 * t51) * t66) * t28 ((-t49 * t51 + t50 * t54) * t30 - (-t49 * t54 - t50 * t51) * t66) * t28 * t63, 0, 0, 0, t57 * t28;];
Ja_rot  = t1;
