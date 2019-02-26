% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:23
% EndTime: 2019-02-26 21:42:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (944->39), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->56)
t88 = cos(pkin(6));
t89 = sin(qJ(2));
t92 = cos(qJ(1));
t96 = t92 * t89;
t90 = sin(qJ(1));
t91 = cos(qJ(2));
t97 = t90 * t91;
t76 = t88 * t96 + t97;
t86 = pkin(11) + qJ(4);
t82 = sin(t86);
t84 = cos(t86);
t87 = sin(pkin(6));
t99 = t87 * t92;
t65 = t76 * t82 + t84 * t99;
t102 = t87 * t89;
t73 = t82 * t102 - t88 * t84;
t64 = atan2(-t65, t73);
t61 = sin(t64);
t62 = cos(t64);
t55 = -t61 * t65 + t62 * t73;
t54 = 0.1e1 / t55 ^ 2;
t101 = t87 * t90;
t95 = t92 * t91;
t98 = t90 * t89;
t78 = -t88 * t98 + t95;
t69 = -t84 * t101 + t78 * t82;
t108 = t54 * t69;
t70 = t82 * t101 + t78 * t84;
t77 = t88 * t97 + t96;
t85 = pkin(12) + qJ(6);
t81 = sin(t85);
t83 = cos(t85);
t60 = t70 * t83 + t77 * t81;
t58 = 0.1e1 / t60 ^ 2;
t59 = t70 * t81 - t77 * t83;
t107 = t58 * t59;
t106 = t62 * t65;
t72 = 0.1e1 / t73 ^ 2;
t105 = t65 * t72;
t104 = t69 ^ 2 * t54;
t103 = t77 * t84;
t100 = t87 * t91;
t94 = t59 ^ 2 * t58 + 0.1e1;
t67 = t76 * t84 - t82 * t99;
t93 = -t61 * t73 - t106;
t75 = t88 * t95 - t98;
t74 = t84 * t102 + t88 * t82;
t71 = 0.1e1 / t73;
t63 = 0.1e1 / (t65 ^ 2 * t72 + 0.1e1);
t57 = 0.1e1 / t60;
t56 = 0.1e1 / t94;
t53 = 0.1e1 / t55;
t52 = 0.1e1 / (0.1e1 + t104);
t51 = (t100 * t105 - t71 * t75) * t82 * t63;
t50 = (t74 * t105 - t67 * t71) * t63;
t1 = [-t69 * t71 * t63, t51, 0, t50, 0, 0; (-t65 * t53 - (-t61 + (t71 * t106 + t61) * t63) * t104) * t52 (-t77 * t82 * t53 - ((t62 * t100 - t61 * t75) * t82 + t93 * t51) * t108) * t52, 0 (t70 * t53 - (t93 * t50 - t61 * t67 + t62 * t74) * t108) * t52, 0, 0; ((-t67 * t81 - t75 * t83) * t57 - (-t67 * t83 + t75 * t81) * t107) * t56 ((-t81 * t103 - t78 * t83) * t57 - (-t83 * t103 + t78 * t81) * t107) * t56, 0 (t83 * t107 - t81 * t57) * t69 * t56, 0, t94 * t56;];
Ja_rot  = t1;
