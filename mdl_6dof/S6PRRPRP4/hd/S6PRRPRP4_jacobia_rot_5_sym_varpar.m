% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:59
% EndTime: 2019-02-26 20:02:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (334->30), mult. (949->78), div. (65->9), fcn. (1349->13), ass. (0->50)
t64 = sin(pkin(10));
t66 = cos(pkin(10));
t73 = cos(qJ(2));
t67 = cos(pkin(6));
t70 = sin(qJ(2));
t77 = t67 * t70;
t58 = t64 * t73 + t66 * t77;
t72 = cos(qJ(3));
t65 = sin(pkin(6));
t69 = sin(qJ(3));
t80 = t65 * t69;
t51 = t58 * t72 - t66 * t80;
t79 = t65 * t72;
t62 = t67 * t69 + t70 * t79;
t49 = atan2(-t51, t62);
t46 = sin(t49);
t47 = cos(t49);
t40 = -t46 * t51 + t47 * t62;
t39 = 0.1e1 / t40 ^ 2;
t60 = -t64 * t77 + t66 * t73;
t54 = t60 * t72 + t64 * t80;
t85 = t39 * t54;
t53 = t60 * t69 - t64 * t79;
t68 = sin(qJ(5));
t76 = t67 * t73;
t59 = t64 * t76 + t66 * t70;
t71 = cos(qJ(5));
t81 = t59 * t71;
t45 = t53 * t68 + t81;
t43 = 0.1e1 / t45 ^ 2;
t82 = t59 * t68;
t44 = -t53 * t71 + t82;
t84 = t43 * t44;
t56 = 0.1e1 / t62 ^ 2;
t83 = t51 * t56;
t78 = t65 * t73;
t75 = t43 * t44 ^ 2 + 0.1e1;
t74 = -t46 * t62 - t47 * t51;
t61 = t67 * t72 - t70 * t80;
t57 = -t64 * t70 + t66 * t76;
t55 = 0.1e1 / t62;
t50 = t58 * t69 + t66 * t79;
t48 = 0.1e1 / (t51 ^ 2 * t56 + 0.1e1);
t42 = 0.1e1 / t45;
t41 = 0.1e1 / t75;
t38 = 0.1e1 / t40;
t37 = 0.1e1 / (t39 * t54 ^ 2 + 0.1e1);
t36 = (-t55 * t57 + t78 * t83) * t72 * t48;
t35 = (t50 * t55 + t61 * t83) * t48;
t1 = [0, t36, t35, 0, 0, 0; 0 (-t59 * t72 * t38 - ((-t46 * t57 + t47 * t78) * t72 + t74 * t36) * t85) * t37 (-t53 * t38 - (t35 * t74 + t46 * t50 + t47 * t61) * t85) * t37, 0, 0, 0; 0 ((t60 * t68 + t69 * t81) * t42 - (t60 * t71 - t69 * t82) * t84) * t41 (-t42 * t71 - t68 * t84) * t54 * t41, 0, t75 * t41, 0;];
Ja_rot  = t1;
