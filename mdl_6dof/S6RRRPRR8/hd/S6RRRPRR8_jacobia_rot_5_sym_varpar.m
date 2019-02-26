% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR8_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:47
% EndTime: 2019-02-26 22:19:48
% DurationCPUTime: 0.23s
% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
t85 = cos(pkin(6));
t88 = sin(qJ(1));
t90 = cos(qJ(2));
t95 = t88 * t90;
t87 = sin(qJ(2));
t91 = cos(qJ(1));
t97 = t87 * t91;
t76 = t85 * t97 + t95;
t83 = qJ(3) + pkin(12);
t81 = sin(t83);
t82 = cos(t83);
t84 = sin(pkin(6));
t98 = t84 * t91;
t65 = t76 * t81 + t82 * t98;
t101 = t84 * t87;
t73 = t81 * t101 - t82 * t85;
t64 = atan2(-t65, t73);
t57 = sin(t64);
t58 = cos(t64);
t55 = -t57 * t65 + t58 * t73;
t54 = 0.1e1 / t55 ^ 2;
t100 = t84 * t88;
t94 = t90 * t91;
t96 = t88 * t87;
t78 = -t85 * t96 + t94;
t69 = -t82 * t100 + t78 * t81;
t108 = t54 * t69;
t107 = t54 * t69 ^ 2;
t106 = t58 * t65;
t77 = t85 * t95 + t97;
t86 = sin(qJ(5));
t103 = t77 * t86;
t70 = t81 * t100 + t78 * t82;
t89 = cos(qJ(5));
t63 = t70 * t89 + t103;
t60 = 0.1e1 / t63 ^ 2;
t102 = t77 * t89;
t62 = t70 * t86 - t102;
t105 = t60 * t62;
t72 = 0.1e1 / t73 ^ 2;
t104 = t65 * t72;
t99 = t84 * t90;
t93 = t60 * t62 ^ 2 + 0.1e1;
t67 = t76 * t82 - t81 * t98;
t92 = -t57 * t73 - t106;
t75 = t85 * t94 - t96;
t74 = t82 * t101 + t81 * t85;
t71 = 0.1e1 / t73;
t61 = 0.1e1 / (t65 ^ 2 * t72 + 0.1e1);
t59 = 0.1e1 / t63;
t56 = 0.1e1 / t93;
t53 = 0.1e1 / t55;
t52 = 0.1e1 / (0.1e1 + t107);
t51 = (t99 * t104 - t71 * t75) * t81 * t61;
t50 = (t74 * t104 - t67 * t71) * t61;
t1 = [-t69 * t71 * t61, t51, t50, 0, 0, 0; (-t65 * t53 - (-t57 + (t71 * t106 + t57) * t61) * t107) * t52 (-t77 * t81 * t53 - ((-t57 * t75 + t58 * t99) * t81 + t92 * t51) * t108) * t52 (t70 * t53 - (t92 * t50 - t57 * t67 + t58 * t74) * t108) * t52, 0, 0, 0; ((-t67 * t86 - t75 * t89) * t59 - (-t67 * t89 + t75 * t86) * t105) * t56 ((-t82 * t103 - t78 * t89) * t59 - (-t82 * t102 + t78 * t86) * t105) * t56 (t89 * t105 - t59 * t86) * t69 * t56, 0, t93 * t56, 0;];
Ja_rot  = t1;
