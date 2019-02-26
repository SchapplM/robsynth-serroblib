% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR3
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR3_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:24
% EndTime: 2019-02-26 20:05:24
% DurationCPUTime: 0.13s
% Computational Cost: add. (239->26), mult. (687->64), div. (36->9), fcn. (960->15), ass. (0->42)
t66 = sin(pkin(12));
t67 = sin(pkin(7));
t68 = sin(pkin(6));
t82 = t67 * t68;
t84 = t66 * t82;
t71 = cos(pkin(7));
t65 = sin(pkin(13));
t69 = cos(pkin(13));
t73 = sin(qJ(3));
t75 = cos(qJ(3));
t77 = t75 * t65 + t73 * t69;
t55 = t77 * t71;
t70 = cos(pkin(12));
t74 = sin(qJ(2));
t72 = cos(pkin(6));
t76 = cos(qJ(2));
t78 = t72 * t76;
t59 = -t66 * t78 - t70 * t74;
t79 = t72 * t74;
t60 = -t66 * t79 + t70 * t76;
t61 = t73 * t65 - t75 * t69;
t46 = t59 * t55 - t60 * t61 + t77 * t84;
t43 = 0.1e1 / t46 ^ 2;
t54 = t61 * t71;
t44 = -t59 * t54 - t60 * t77 - t61 * t84;
t83 = t44 ^ 2 * t43;
t81 = t68 * t71;
t80 = t68 * t74;
t58 = -t66 * t76 - t70 * t79;
t57 = t72 * t71 - t76 * t82;
t56 = 0.1e1 / t57 ^ 2;
t52 = -t59 * t67 + t66 * t81;
t51 = (-t66 * t74 + t70 * t78) * t67 + t70 * t81;
t50 = atan2(t51, t57);
t48 = cos(t50);
t47 = sin(t50);
t42 = 0.1e1 / t46;
t41 = t47 * t51 + t48 * t57;
t40 = 0.1e1 / t41 ^ 2;
t38 = 0.1e1 / (0.1e1 + t83);
t37 = (t58 / t57 - t51 * t56 * t80) * t67 / (t51 ^ 2 * t56 + 0.1e1);
t1 = [0, t37, 0, 0, 0, 0; 0 (t60 * t67 / t41 - ((t47 * t58 + t48 * t80) * t67 + (-t47 * t57 + t48 * t51) * t37) * t52 * t40) / (t52 ^ 2 * t40 + 0.1e1) 0, 0, 0, 0; 0 ((-t60 * t54 + t59 * t77) * t42 + (-t60 * t55 - t59 * t61) * t44 * t43) * t38 (t46 * t42 + t83) * t38, 0, 0, 0;];
Ja_rot  = t1;
