% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:11
% EndTime: 2019-02-26 20:01:11
% DurationCPUTime: 0.12s
% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
t79 = sin(pkin(10));
t81 = cos(pkin(10));
t86 = cos(qJ(2));
t82 = cos(pkin(6));
t84 = sin(qJ(2));
t90 = t82 * t84;
t72 = t79 * t86 + t81 * t90;
t78 = qJ(3) + pkin(11);
t76 = sin(t78);
t77 = cos(t78);
t80 = sin(pkin(6));
t93 = t80 * t81;
t62 = t72 * t76 + t77 * t93;
t92 = t80 * t84;
t69 = t76 * t92 - t82 * t77;
t61 = atan2(-t62, t69);
t56 = sin(t61);
t57 = cos(t61);
t52 = -t56 * t62 + t57 * t69;
t51 = 0.1e1 / t52 ^ 2;
t74 = -t79 * t90 + t81 * t86;
t94 = t79 * t80;
t65 = t74 * t76 - t77 * t94;
t99 = t51 * t65;
t66 = t74 * t77 + t76 * t94;
t85 = cos(qJ(5));
t89 = t82 * t86;
t73 = t79 * t89 + t81 * t84;
t83 = sin(qJ(5));
t96 = t73 * t83;
t59 = t66 * t85 + t96;
t55 = 0.1e1 / t59 ^ 2;
t95 = t73 * t85;
t58 = t66 * t83 - t95;
t98 = t55 * t58;
t68 = 0.1e1 / t69 ^ 2;
t97 = t62 * t68;
t91 = t80 * t86;
t88 = t58 ^ 2 * t55 + 0.1e1;
t87 = -t56 * t69 - t57 * t62;
t71 = -t79 * t84 + t81 * t89;
t70 = t82 * t76 + t77 * t92;
t67 = 0.1e1 / t69;
t64 = t72 * t77 - t76 * t93;
t60 = 0.1e1 / (t62 ^ 2 * t68 + 0.1e1);
t54 = 0.1e1 / t59;
t53 = 0.1e1 / t88;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (t65 ^ 2 * t51 + 0.1e1);
t48 = (-t67 * t71 + t91 * t97) * t76 * t60;
t47 = (-t64 * t67 + t70 * t97) * t60;
t1 = [0, t48, t47, 0, 0, 0; 0 (-t73 * t76 * t50 - ((-t56 * t71 + t57 * t91) * t76 + t87 * t48) * t99) * t49 (t66 * t50 - (t87 * t47 - t56 * t64 + t57 * t70) * t99) * t49, 0, 0, 0; 0 ((-t74 * t85 - t77 * t96) * t54 - (t74 * t83 - t77 * t95) * t98) * t53 (-t54 * t83 + t85 * t98) * t65 * t53, 0, t88 * t53, 0;];
Ja_rot  = t1;
