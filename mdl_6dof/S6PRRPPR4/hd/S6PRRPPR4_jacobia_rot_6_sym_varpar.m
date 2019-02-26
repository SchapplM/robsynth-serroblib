% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR4
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
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:00
% EndTime: 2019-02-26 20:00:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (453->35), mult. (1262->91), div. (65->9), fcn. (1781->15), ass. (0->55)
t88 = sin(pkin(6));
t96 = cos(qJ(3));
t103 = t88 * t96;
t91 = cos(pkin(6));
t94 = sin(qJ(2));
t101 = t91 * t94;
t87 = sin(pkin(10));
t90 = cos(pkin(10));
t97 = cos(qJ(2));
t81 = t90 * t101 + t87 * t97;
t93 = sin(qJ(3));
t74 = t90 * t103 + t81 * t93;
t104 = t88 * t93;
t84 = -t94 * t104 + t91 * t96;
t71 = atan2(t74, t84);
t68 = sin(t71);
t69 = cos(t71);
t61 = t68 * t74 + t69 * t84;
t60 = 0.1e1 / t61 ^ 2;
t83 = -t87 * t101 + t90 * t97;
t76 = t87 * t103 - t83 * t93;
t108 = t60 * t76;
t77 = t87 * t104 + t83 * t96;
t100 = t91 * t97;
t82 = t87 * t100 + t90 * t94;
t86 = sin(pkin(11));
t89 = cos(pkin(11));
t66 = t77 * t86 - t82 * t89;
t67 = t77 * t89 + t82 * t86;
t92 = sin(qJ(6));
t95 = cos(qJ(6));
t65 = t66 * t92 + t67 * t95;
t63 = 0.1e1 / t65 ^ 2;
t64 = -t66 * t95 + t67 * t92;
t107 = t63 * t64;
t79 = 0.1e1 / t84 ^ 2;
t106 = t74 * t79;
t105 = t82 * t96;
t102 = t88 * t97;
t99 = t64 ^ 2 * t63 + 0.1e1;
t98 = -t68 * t84 + t69 * t74;
t85 = -t94 * t103 - t91 * t93;
t80 = t90 * t100 - t87 * t94;
t78 = 0.1e1 / t84;
t75 = -t90 * t104 + t81 * t96;
t73 = -t89 * t105 + t83 * t86;
t72 = -t86 * t105 - t83 * t89;
t70 = 0.1e1 / (t74 ^ 2 * t79 + 0.1e1);
t62 = 0.1e1 / t65;
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (t76 ^ 2 * t60 + 0.1e1);
t57 = (t102 * t106 + t78 * t80) * t93 * t70;
t56 = (-t85 * t106 + t75 * t78) * t70;
t55 = 0.1e1 / t99;
t1 = [0, t57, t56, 0, 0, 0; 0 (t82 * t93 * t59 - ((-t69 * t102 + t68 * t80) * t93 + t98 * t57) * t108) * t58 (-t77 * t59 - (t98 * t56 + t68 * t75 + t69 * t85) * t108) * t58, 0, 0, 0; 0 ((-t72 * t95 + t73 * t92) * t62 - (t72 * t92 + t73 * t95) * t107) * t55 ((-t86 * t95 + t89 * t92) * t62 - (t86 * t92 + t89 * t95) * t107) * t55 * t76, 0, 0, t99 * t55;];
Ja_rot  = t1;
