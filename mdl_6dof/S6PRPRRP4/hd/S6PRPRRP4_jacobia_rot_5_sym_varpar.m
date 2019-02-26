% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRP4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:01
% EndTime: 2019-02-26 19:52:01
% DurationCPUTime: 0.14s
% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
t68 = sin(pkin(10));
t70 = cos(pkin(10));
t75 = cos(qJ(2));
t71 = cos(pkin(6));
t73 = sin(qJ(2));
t79 = t71 * t73;
t61 = t68 * t75 + t70 * t79;
t67 = pkin(11) + qJ(4);
t65 = sin(t67);
t66 = cos(t67);
t69 = sin(pkin(6));
t82 = t69 * t70;
t51 = t61 * t65 + t66 * t82;
t81 = t69 * t73;
t58 = t65 * t81 - t66 * t71;
t50 = atan2(-t51, t58);
t45 = sin(t50);
t46 = cos(t50);
t41 = -t45 * t51 + t46 * t58;
t40 = 0.1e1 / t41 ^ 2;
t63 = -t68 * t79 + t70 * t75;
t83 = t68 * t69;
t54 = t63 * t65 - t66 * t83;
t88 = t40 * t54;
t55 = t63 * t66 + t65 * t83;
t74 = cos(qJ(5));
t78 = t71 * t75;
t62 = t68 * t78 + t70 * t73;
t72 = sin(qJ(5));
t85 = t62 * t72;
t48 = t55 * t74 + t85;
t44 = 0.1e1 / t48 ^ 2;
t84 = t62 * t74;
t47 = t55 * t72 - t84;
t87 = t44 * t47;
t57 = 0.1e1 / t58 ^ 2;
t86 = t51 * t57;
t80 = t69 * t75;
t77 = t44 * t47 ^ 2 + 0.1e1;
t76 = -t45 * t58 - t46 * t51;
t60 = -t68 * t73 + t70 * t78;
t59 = t65 * t71 + t66 * t81;
t56 = 0.1e1 / t58;
t53 = t61 * t66 - t65 * t82;
t49 = 0.1e1 / (t51 ^ 2 * t57 + 0.1e1);
t43 = 0.1e1 / t48;
t42 = 0.1e1 / t77;
t39 = 0.1e1 / t41;
t38 = 0.1e1 / (t40 * t54 ^ 2 + 0.1e1);
t37 = (-t56 * t60 + t80 * t86) * t65 * t49;
t36 = (-t53 * t56 + t59 * t86) * t49;
t1 = [0, t37, 0, t36, 0, 0; 0 (-t62 * t65 * t39 - ((-t45 * t60 + t46 * t80) * t65 + t76 * t37) * t88) * t38, 0 (t55 * t39 - (t36 * t76 - t45 * t53 + t46 * t59) * t88) * t38, 0, 0; 0 ((-t63 * t74 - t66 * t85) * t43 - (t63 * t72 - t66 * t84) * t87) * t42, 0 (-t43 * t72 + t74 * t87) * t54 * t42, t77 * t42, 0;];
Ja_rot  = t1;
