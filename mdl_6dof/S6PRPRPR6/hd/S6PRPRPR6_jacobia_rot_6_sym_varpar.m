% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:12
% EndTime: 2019-02-26 19:49:12
% DurationCPUTime: 0.13s
% Computational Cost: add. (382->30), mult. (949->80), div. (65->9), fcn. (1349->13), ass. (0->51)
t69 = sin(pkin(10));
t71 = cos(pkin(10));
t74 = sin(qJ(2));
t72 = cos(pkin(6));
t76 = cos(qJ(2));
t79 = t72 * t76;
t60 = t69 * t74 - t71 * t79;
t75 = cos(qJ(4));
t70 = sin(pkin(6));
t73 = sin(qJ(4));
t84 = t70 * t73;
t55 = t60 * t75 + t71 * t84;
t81 = t70 * t76;
t64 = t72 * t73 + t75 * t81;
t52 = atan2(t55, t64);
t49 = sin(t52);
t50 = cos(t52);
t43 = t49 * t55 + t50 * t64;
t42 = 0.1e1 / t43 ^ 2;
t62 = t69 * t79 + t71 * t74;
t53 = -t62 * t75 + t69 * t84;
t88 = t42 * t53;
t82 = t70 * t75;
t54 = t62 * t73 + t69 * t82;
t80 = t72 * t74;
t63 = -t69 * t80 + t71 * t76;
t68 = pkin(11) + qJ(6);
t66 = sin(t68);
t67 = cos(t68);
t48 = t54 * t67 + t63 * t66;
t46 = 0.1e1 / t48 ^ 2;
t47 = t54 * t66 - t63 * t67;
t87 = t46 * t47;
t59 = 0.1e1 / t64 ^ 2;
t86 = t55 * t59;
t85 = t63 * t73;
t83 = t70 * t74;
t78 = t46 * t47 ^ 2 + 0.1e1;
t77 = -t49 * t64 + t50 * t55;
t65 = t72 * t75 - t73 * t81;
t61 = t69 * t76 + t71 * t80;
t58 = 0.1e1 / t64;
t56 = -t60 * t73 + t71 * t82;
t51 = 0.1e1 / (t55 ^ 2 * t59 + 0.1e1);
t45 = 0.1e1 / t48;
t44 = 0.1e1 / t78;
t41 = 0.1e1 / t43;
t40 = 0.1e1 / (t42 * t53 ^ 2 + 0.1e1);
t39 = (t58 * t61 + t83 * t86) * t75 * t51;
t38 = (t56 * t58 - t65 * t86) * t51;
t1 = [0, t39, 0, t38, 0, 0; 0 (-t63 * t75 * t41 - ((t49 * t61 - t50 * t83) * t75 + t77 * t39) * t88) * t40, 0 (t54 * t41 - (t38 * t77 + t49 * t56 + t50 * t65) * t88) * t40, 0, 0; 0 ((t62 * t67 + t66 * t85) * t45 - (-t62 * t66 + t67 * t85) * t87) * t44, 0 (-t45 * t66 + t67 * t87) * t53 * t44, 0, t78 * t44;];
Ja_rot  = t1;
