% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:01
% EndTime: 2019-02-26 22:07:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (209->26), mult. (425->68), div. (58->9), fcn. (611->11), ass. (0->41)
t67 = sin(qJ(1));
t68 = cos(qJ(3));
t69 = cos(qJ(2));
t65 = sin(qJ(3));
t70 = cos(qJ(1));
t73 = t70 * t65;
t52 = -t67 * t68 + t69 * t73;
t72 = t70 * t68;
t53 = t67 * t65 + t69 * t72;
t61 = pkin(10) + qJ(6);
t59 = sin(t61);
t60 = cos(t61);
t44 = t52 * t59 + t53 * t60;
t42 = 0.1e1 / t44 ^ 2;
t43 = -t52 * t60 + t53 * t59;
t81 = t42 * t43;
t80 = t43 ^ 2 * t42;
t66 = sin(qJ(2));
t75 = t67 * t66;
t57 = atan2(t75, t69);
t54 = sin(t57);
t55 = cos(t57);
t48 = t54 * t75 + t55 * t69;
t47 = 0.1e1 / t48 ^ 2;
t79 = t47 * t70 ^ 2;
t62 = t66 ^ 2;
t78 = t62 / t69 ^ 2;
t77 = t66 * t70;
t56 = 0.1e1 / (t67 ^ 2 * t78 + 0.1e1);
t76 = t67 * t56;
t74 = t67 * t69;
t71 = 0.1e1 + t80;
t63 = 0.1e1 / t69;
t51 = -t68 * t74 + t73;
t50 = -t65 * t74 - t72;
t49 = (0.1e1 + t78) * t76;
t46 = 0.1e1 / t48;
t45 = 0.1e1 / (t62 * t79 + 0.1e1);
t41 = 0.1e1 / t44;
t40 = 0.1e1 / t71;
t1 = [t63 * t56 * t77, t49, 0, 0, 0, 0; (t46 * t75 + (t55 * t62 * t63 * t76 + (-t56 + 0.1e1) * t66 * t54) * t66 * t79) * t45 (-t69 * t46 + (t54 * t74 - t55 * t66 + (-t54 * t69 + t55 * t75) * t49) * t66 * t47) * t70 * t45, 0, 0, 0, 0; ((-t50 * t60 + t51 * t59) * t41 - (t50 * t59 + t51 * t60) * t81) * t40 ((-t59 * t68 + t60 * t65) * t41 - (-t59 * t65 - t60 * t68) * t81) * t40 * t77 (-t44 * t41 - t80) * t40, 0, 0, t71 * t40;];
Ja_rot  = t1;
