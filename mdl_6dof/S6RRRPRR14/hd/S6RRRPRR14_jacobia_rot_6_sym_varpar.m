% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR14_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:40
% EndTime: 2019-02-26 22:23:40
% DurationCPUTime: 0.22s
% Computational Cost: add. (570->38), mult. (1383->91), div. (90->9), fcn. (1970->13), ass. (0->56)
t90 = sin(pkin(6));
t97 = cos(qJ(1));
t104 = t90 * t97;
t93 = sin(qJ(2));
t101 = t97 * t93;
t94 = sin(qJ(1));
t96 = cos(qJ(2));
t102 = t94 * t96;
t91 = cos(pkin(6));
t83 = t91 * t101 + t102;
t92 = sin(qJ(3));
t95 = cos(qJ(3));
t73 = -t92 * t104 + t83 * t95;
t106 = t90 * t95;
t81 = t93 * t106 + t91 * t92;
t71 = atan2(-t73, t81);
t68 = sin(t71);
t69 = cos(t71);
t62 = -t68 * t73 + t69 * t81;
t61 = 0.1e1 / t62 ^ 2;
t107 = t90 * t92;
t100 = t97 * t96;
t103 = t94 * t93;
t85 = -t91 * t103 + t100;
t77 = t94 * t107 + t85 * t95;
t113 = t61 * t77;
t76 = -t94 * t106 + t85 * t92;
t84 = t91 * t102 + t101;
t89 = qJ(5) + qJ(6);
t87 = sin(t89);
t88 = cos(t89);
t67 = t76 * t87 + t84 * t88;
t65 = 0.1e1 / t67 ^ 2;
t66 = -t76 * t88 + t84 * t87;
t112 = t65 * t66;
t111 = t69 * t73;
t79 = 0.1e1 / t81 ^ 2;
t110 = t73 * t79;
t109 = t77 ^ 2 * t61;
t108 = t84 * t92;
t105 = t90 * t96;
t99 = t66 ^ 2 * t65 + 0.1e1;
t98 = -t68 * t81 - t111;
t72 = t95 * t104 + t83 * t92;
t82 = t91 * t100 - t103;
t80 = -t93 * t107 + t91 * t95;
t78 = 0.1e1 / t81;
t70 = 0.1e1 / (t73 ^ 2 * t79 + 0.1e1);
t64 = 0.1e1 / t67;
t63 = 0.1e1 / t99;
t60 = 0.1e1 / t62;
t59 = 0.1e1 / (0.1e1 + t109);
t58 = (t105 * t110 - t78 * t82) * t95 * t70;
t57 = (t80 * t110 + t72 * t78) * t70;
t56 = t99 * t63;
t1 = [-t77 * t78 * t70, t58, t57, 0, 0, 0; (-t73 * t60 - (-t68 + (t78 * t111 + t68) * t70) * t109) * t59 (-t84 * t95 * t60 - ((t69 * t105 - t68 * t82) * t95 + t98 * t58) * t113) * t59 (-t76 * t60 - (t98 * t57 + t68 * t72 + t69 * t80) * t113) * t59, 0, 0, 0; ((t72 * t88 + t82 * t87) * t64 - (-t72 * t87 + t82 * t88) * t112) * t63 ((t88 * t108 + t85 * t87) * t64 - (-t87 * t108 + t85 * t88) * t112) * t63 (-t87 * t112 - t88 * t64) * t77 * t63, 0, t56, t56;];
Ja_rot  = t1;
