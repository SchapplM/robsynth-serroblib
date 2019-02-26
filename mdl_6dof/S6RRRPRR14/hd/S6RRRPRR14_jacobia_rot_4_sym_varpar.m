% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR14_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:40
% EndTime: 2019-02-26 22:23:40
% DurationCPUTime: 0.15s
% Computational Cost: add. (355->31), mult. (1022->76), div. (77->9), fcn. (1479->11), ass. (0->48)
t63 = cos(pkin(6));
t65 = sin(qJ(2));
t69 = cos(qJ(1));
t72 = t69 * t65;
t66 = sin(qJ(1));
t68 = cos(qJ(2));
t73 = t66 * t68;
t57 = t63 * t72 + t73;
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t62 = sin(pkin(6));
t75 = t62 * t69;
t45 = t57 * t64 + t67 * t75;
t78 = t62 * t64;
t54 = -t63 * t67 + t65 * t78;
t44 = atan2(-t45, t54);
t40 = sin(t44);
t41 = cos(t44);
t39 = -t40 * t45 + t41 * t54;
t38 = 0.1e1 / t39 ^ 2;
t71 = t69 * t68;
t74 = t66 * t65;
t59 = -t63 * t74 + t71;
t77 = t62 * t67;
t48 = t59 * t64 - t66 * t77;
t83 = t38 * t48;
t82 = t41 * t45;
t51 = 0.1e1 / t54 ^ 2;
t81 = t45 * t51;
t80 = t48 ^ 2 * t38;
t49 = t59 * t67 + t66 * t78;
t58 = t63 * t73 + t72;
t53 = 0.1e1 / t58 ^ 2;
t79 = t49 * t53;
t76 = t62 * t68;
t47 = t57 * t67 - t64 * t75;
t70 = -t40 * t54 - t82;
t56 = t63 * t71 - t74;
t55 = t63 * t64 + t65 * t77;
t52 = 0.1e1 / t58;
t50 = 0.1e1 / t54;
t43 = 0.1e1 / (t49 ^ 2 * t53 + 0.1e1);
t42 = 0.1e1 / (t45 ^ 2 * t51 + 0.1e1);
t37 = 0.1e1 / t39;
t36 = 0.1e1 / (0.1e1 + t80);
t35 = (-t50 * t56 + t76 * t81) * t64 * t42;
t34 = (-t47 * t50 + t55 * t81) * t42;
t1 = [-t48 * t50 * t42, t35, t34, 0, 0, 0; (-t45 * t37 - (-t40 + (t50 * t82 + t40) * t42) * t80) * t36 (-t58 * t64 * t37 - ((-t40 * t56 + t41 * t76) * t64 + t70 * t35) * t83) * t36 (t49 * t37 - (t70 * t34 - t40 * t47 + t41 * t55) * t83) * t36, 0, 0, 0; (-t47 * t52 - t56 * t79) * t43 (-t52 * t58 * t67 - t59 * t79) * t43, -t48 * t52 * t43, 0, 0, 0;];
Ja_rot  = t1;
