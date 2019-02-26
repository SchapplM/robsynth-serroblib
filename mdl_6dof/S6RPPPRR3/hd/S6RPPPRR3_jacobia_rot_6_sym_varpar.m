% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPPRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:23:46
% EndTime: 2019-02-26 20:23:47
% DurationCPUTime: 0.09s
% Computational Cost: add. (325->21), mult. (442->61), div. (52->9), fcn. (659->11), ass. (0->39)
t59 = sin(pkin(9));
t60 = cos(pkin(9));
t70 = sin(qJ(1));
t71 = cos(qJ(1));
t43 = -t70 * t59 - t71 * t60;
t72 = t43 ^ 2;
t55 = pkin(10) + qJ(5);
t54 = cos(t55);
t44 = t71 * t59 - t70 * t60;
t53 = sin(t55);
t64 = t44 * t53;
t41 = atan2(t64, t54);
t39 = sin(t41);
t40 = cos(t41);
t33 = t39 * t64 + t40 * t54;
t32 = 0.1e1 / t33 ^ 2;
t69 = t32 * t53;
t56 = sin(qJ(6));
t57 = cos(qJ(6));
t61 = t54 * t57;
t38 = -t43 * t61 + t44 * t56;
t36 = 0.1e1 / t38 ^ 2;
t62 = t54 * t56;
t37 = -t43 * t62 - t44 * t57;
t68 = t36 * t37;
t67 = t39 * t54;
t66 = t43 * t53;
t50 = t53 ^ 2;
t63 = t50 / t54 ^ 2;
t42 = 0.1e1 / (t44 ^ 2 * t63 + 0.1e1);
t65 = t44 * t42;
t58 = t37 ^ 2 * t36 + 0.1e1;
t51 = 0.1e1 / t54;
t35 = 0.1e1 / t38;
t34 = 0.1e1 / t58;
t31 = 0.1e1 / t33;
t30 = (0.1e1 + t63) * t65;
t29 = 0.1e1 / (t72 * t50 * t32 + 0.1e1);
t1 = [t51 * t42 * t66, 0, 0, 0, t30, 0; (t31 * t64 + (t40 * t50 * t51 * t65 + (-t42 + 0.1e1) * t53 * t39) * t72 * t69) * t29, 0, 0, 0 (-t54 * t31 + (t44 * t67 - t40 * t53 + (t40 * t64 - t67) * t30) * t69) * t43 * t29, 0; ((-t43 * t57 + t44 * t62) * t35 - (t43 * t56 + t44 * t61) * t68) * t34, 0, 0, 0 (t35 * t56 - t57 * t68) * t34 * t66, t58 * t34;];
Ja_rot  = t1;
