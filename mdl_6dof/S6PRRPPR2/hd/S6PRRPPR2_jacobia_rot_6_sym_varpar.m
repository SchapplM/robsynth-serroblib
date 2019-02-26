% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:46
% EndTime: 2019-02-26 19:58:46
% DurationCPUTime: 0.12s
% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
t72 = sin(pkin(10));
t74 = cos(pkin(10));
t79 = cos(qJ(2));
t75 = cos(pkin(6));
t77 = sin(qJ(2));
t83 = t75 * t77;
t65 = t72 * t79 + t74 * t83;
t71 = qJ(3) + pkin(11);
t69 = sin(t71);
t70 = cos(t71);
t73 = sin(pkin(6));
t86 = t73 * t74;
t56 = t65 * t70 - t69 * t86;
t85 = t73 * t77;
t63 = t75 * t69 + t70 * t85;
t54 = atan2(-t56, t63);
t49 = sin(t54);
t50 = cos(t54);
t45 = -t49 * t56 + t50 * t63;
t44 = 0.1e1 / t45 ^ 2;
t67 = -t72 * t83 + t74 * t79;
t87 = t72 * t73;
t59 = t67 * t70 + t69 * t87;
t92 = t44 * t59;
t58 = t67 * t69 - t70 * t87;
t76 = sin(qJ(6));
t82 = t75 * t79;
t66 = t72 * t82 + t74 * t77;
t78 = cos(qJ(6));
t88 = t66 * t78;
t52 = t58 * t76 + t88;
t48 = 0.1e1 / t52 ^ 2;
t89 = t66 * t76;
t51 = -t58 * t78 + t89;
t91 = t48 * t51;
t61 = 0.1e1 / t63 ^ 2;
t90 = t56 * t61;
t84 = t73 * t79;
t81 = t51 ^ 2 * t48 + 0.1e1;
t80 = -t49 * t63 - t50 * t56;
t64 = -t72 * t77 + t74 * t82;
t62 = -t69 * t85 + t75 * t70;
t60 = 0.1e1 / t63;
t55 = t65 * t69 + t70 * t86;
t53 = 0.1e1 / (t56 ^ 2 * t61 + 0.1e1);
t47 = 0.1e1 / t52;
t46 = 0.1e1 / t81;
t43 = 0.1e1 / t45;
t42 = 0.1e1 / (t59 ^ 2 * t44 + 0.1e1);
t41 = (-t60 * t64 + t84 * t90) * t70 * t53;
t40 = (t55 * t60 + t62 * t90) * t53;
t1 = [0, t41, t40, 0, 0, 0; 0 (-t66 * t70 * t43 - ((-t49 * t64 + t50 * t84) * t70 + t80 * t41) * t92) * t42 (-t58 * t43 - (t80 * t40 + t49 * t55 + t50 * t62) * t92) * t42, 0, 0, 0; 0 ((t67 * t76 + t69 * t88) * t47 - (t67 * t78 - t69 * t89) * t91) * t46 (-t47 * t78 - t76 * t91) * t59 * t46, 0, 0, t81 * t46;];
Ja_rot  = t1;
