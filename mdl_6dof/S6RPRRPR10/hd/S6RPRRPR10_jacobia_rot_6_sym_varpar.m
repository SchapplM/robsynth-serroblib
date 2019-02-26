% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:03
% EndTime: 2019-02-26 21:06:03
% DurationCPUTime: 0.10s
% Computational Cost: add. (164->26), mult. (425->68), div. (58->9), fcn. (611->11), ass. (0->43)
t63 = sin(qJ(1));
t81 = t63 ^ 2;
t62 = sin(qJ(3));
t65 = cos(qJ(4));
t67 = cos(qJ(1));
t70 = t67 * t65;
t61 = sin(qJ(4));
t75 = t63 * t61;
t48 = t62 * t75 - t70;
t71 = t67 * t61;
t74 = t63 * t65;
t49 = t62 * t74 + t71;
t60 = sin(qJ(6));
t64 = cos(qJ(6));
t43 = t48 * t60 + t49 * t64;
t41 = 0.1e1 / t43 ^ 2;
t42 = -t48 * t64 + t49 * t60;
t80 = t41 * t42;
t79 = t42 ^ 2 * t41;
t66 = cos(qJ(3));
t69 = t67 * t66;
t54 = atan2(t69, -t62);
t52 = sin(t54);
t53 = cos(t54);
t46 = t52 * t69 - t53 * t62;
t45 = 0.1e1 / t46 ^ 2;
t78 = t45 * t66;
t77 = t52 * t62;
t59 = t66 ^ 2;
t76 = 0.1e1 / t62 ^ 2 * t59;
t73 = t63 * t66;
t55 = 0.1e1 / (t67 ^ 2 * t76 + 0.1e1);
t72 = t67 * t55;
t68 = 0.1e1 + t79;
t57 = 0.1e1 / t62;
t51 = t62 * t70 - t75;
t50 = t62 * t71 + t74;
t47 = (0.1e1 + t76) * t72;
t44 = 0.1e1 / t46;
t40 = 0.1e1 / t43;
t39 = 0.1e1 / (t81 * t59 * t45 + 0.1e1);
t38 = 0.1e1 / t68;
t1 = [t57 * t55 * t73, 0, t47, 0, 0, 0; (t44 * t69 - (t53 * t57 * t59 * t72 + (t55 - 0.1e1) * t66 * t52) * t81 * t78) * t39, 0 (-t62 * t44 - (-t67 * t77 - t53 * t66 + (t53 * t69 + t77) * t47) * t78) * t63 * t39, 0, 0, 0; ((-t50 * t64 + t51 * t60) * t40 - (t50 * t60 + t51 * t64) * t80) * t38, 0 ((t60 * t65 - t61 * t64) * t40 - (t60 * t61 + t64 * t65) * t80) * t38 * t73 (-t43 * t40 - t79) * t38, 0, t68 * t38;];
Ja_rot  = t1;
