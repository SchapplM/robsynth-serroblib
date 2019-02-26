% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRP3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP3_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:25
% EndTime: 2019-02-26 20:16:25
% DurationCPUTime: 0.13s
% Computational Cost: add. (334->30), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->49)
t62 = sin(pkin(11));
t64 = cos(pkin(11));
t71 = cos(qJ(2));
t65 = cos(pkin(6));
t68 = sin(qJ(2));
t75 = t65 * t68;
t56 = t62 * t71 + t64 * t75;
t67 = sin(qJ(3));
t63 = sin(pkin(6));
t70 = cos(qJ(3));
t77 = t63 * t70;
t48 = t56 * t67 + t64 * t77;
t78 = t63 * t67;
t59 = -t65 * t70 + t68 * t78;
t47 = atan2(-t48, t59);
t44 = sin(t47);
t45 = cos(t47);
t38 = -t44 * t48 + t45 * t59;
t37 = 0.1e1 / t38 ^ 2;
t58 = -t62 * t75 + t64 * t71;
t51 = t58 * t67 - t62 * t77;
t82 = t37 * t51;
t52 = t58 * t70 + t62 * t78;
t74 = t65 * t71;
t57 = t62 * t74 + t64 * t68;
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t43 = t52 * t69 + t57 * t66;
t41 = 0.1e1 / t43 ^ 2;
t42 = t52 * t66 - t57 * t69;
t81 = t41 * t42;
t54 = 0.1e1 / t59 ^ 2;
t80 = t48 * t54;
t79 = t57 * t70;
t76 = t63 * t71;
t73 = t42 ^ 2 * t41 + 0.1e1;
t72 = -t44 * t59 - t45 * t48;
t60 = t65 * t67 + t68 * t77;
t55 = -t62 * t68 + t64 * t74;
t53 = 0.1e1 / t59;
t50 = t56 * t70 - t64 * t78;
t46 = 0.1e1 / (t48 ^ 2 * t54 + 0.1e1);
t40 = 0.1e1 / t43;
t39 = 0.1e1 / t73;
t36 = 0.1e1 / t38;
t35 = 0.1e1 / (t51 ^ 2 * t37 + 0.1e1);
t34 = (-t53 * t55 + t76 * t80) * t67 * t46;
t33 = (-t50 * t53 + t60 * t80) * t46;
t1 = [0, t34, t33, 0, 0, 0; 0 (-t57 * t67 * t36 - ((-t44 * t55 + t45 * t76) * t67 + t72 * t34) * t82) * t35 (t52 * t36 - (t72 * t33 - t44 * t50 + t45 * t60) * t82) * t35, 0, 0, 0; 0 ((-t58 * t69 - t66 * t79) * t40 - (t58 * t66 - t69 * t79) * t81) * t39 (-t40 * t66 + t69 * t81) * t51 * t39, t73 * t39, 0, 0;];
Ja_rot  = t1;
