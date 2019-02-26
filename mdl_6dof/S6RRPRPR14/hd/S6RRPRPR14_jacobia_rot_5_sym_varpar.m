% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR14_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:25
% EndTime: 2019-02-26 21:45:26
% DurationCPUTime: 0.13s
% Computational Cost: add. (353->30), mult. (1015->74), div. (76->9), fcn. (1469->11), ass. (0->48)
t64 = cos(pkin(6));
t69 = cos(qJ(2));
t70 = cos(qJ(1));
t73 = t70 * t69;
t66 = sin(qJ(2));
t67 = sin(qJ(1));
t76 = t67 * t66;
t59 = -t64 * t73 + t76;
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t63 = sin(pkin(6));
t77 = t63 * t70;
t50 = t59 * t68 + t65 * t77;
t78 = t63 * t69;
t57 = t64 * t65 + t68 * t78;
t47 = atan2(t50, t57);
t43 = sin(t47);
t44 = cos(t47);
t42 = t43 * t50 + t44 * t57;
t41 = 0.1e1 / t42 ^ 2;
t74 = t70 * t66;
t75 = t67 * t69;
t71 = t64 * t75 + t74;
t79 = t63 * t67;
t48 = t65 * t79 - t71 * t68;
t85 = t41 * t48;
t84 = t44 * t50;
t83 = t48 ^ 2 * t41;
t49 = t71 * t65 + t68 * t79;
t61 = -t64 * t76 + t73;
t56 = 0.1e1 / t61 ^ 2;
t82 = t49 * t56;
t54 = 0.1e1 / t57 ^ 2;
t81 = t50 * t54;
t80 = t63 * t66;
t72 = -t43 * t57 + t84;
t51 = -t59 * t65 + t68 * t77;
t60 = t64 * t74 + t75;
t58 = t64 * t68 - t65 * t78;
t55 = 0.1e1 / t61;
t53 = 0.1e1 / t57;
t46 = 0.1e1 / (t49 ^ 2 * t56 + 0.1e1);
t45 = 0.1e1 / (t50 ^ 2 * t54 + 0.1e1);
t40 = 0.1e1 / t42;
t39 = 0.1e1 / (0.1e1 + t83);
t38 = (t53 * t60 + t80 * t81) * t68 * t45;
t37 = (t51 * t53 - t58 * t81) * t45;
t1 = [-t48 * t53 * t45, t38, 0, t37, 0, 0; (t50 * t40 - (-t43 + (-t53 * t84 + t43) * t45) * t83) * t39 (-t61 * t68 * t40 - ((t43 * t60 - t44 * t80) * t68 + t72 * t38) * t85) * t39, 0 (t49 * t40 - (t72 * t37 + t43 * t51 + t44 * t58) * t85) * t39, 0, 0; (t51 * t55 + t60 * t82) * t46 (t71 * t82 + t65) * t46, 0, -t48 * t55 * t46, 0, 0;];
Ja_rot  = t1;
