% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:44:21
% EndTime: 2019-02-26 22:44:21
% DurationCPUTime: 0.23s
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
t82 = t91 * t101 + t102;
t92 = sin(qJ(3));
t95 = cos(qJ(3));
t71 = t95 * t104 + t82 * t92;
t107 = t90 * t92;
t79 = t93 * t107 - t91 * t95;
t70 = atan2(-t71, t79);
t67 = sin(t70);
t68 = cos(t70);
t61 = -t67 * t71 + t68 * t79;
t60 = 0.1e1 / t61 ^ 2;
t106 = t90 * t95;
t100 = t97 * t96;
t103 = t94 * t93;
t84 = -t91 * t103 + t100;
t75 = -t94 * t106 + t84 * t92;
t113 = t60 * t75;
t76 = t94 * t107 + t84 * t95;
t83 = t91 * t102 + t101;
t89 = qJ(4) + qJ(5);
t87 = sin(t89);
t88 = cos(t89);
t66 = t76 * t88 + t83 * t87;
t64 = 0.1e1 / t66 ^ 2;
t65 = t76 * t87 - t83 * t88;
t112 = t64 * t65;
t111 = t68 * t71;
t78 = 0.1e1 / t79 ^ 2;
t110 = t71 * t78;
t109 = t75 ^ 2 * t60;
t108 = t83 * t95;
t105 = t90 * t96;
t99 = t65 ^ 2 * t64 + 0.1e1;
t73 = -t92 * t104 + t82 * t95;
t98 = -t67 * t79 - t111;
t81 = t91 * t100 - t103;
t80 = t93 * t106 + t91 * t92;
t77 = 0.1e1 / t79;
t69 = 0.1e1 / (t71 ^ 2 * t78 + 0.1e1);
t63 = 0.1e1 / t66;
t62 = 0.1e1 / t99;
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (0.1e1 + t109);
t57 = (t105 * t110 - t77 * t81) * t92 * t69;
t56 = (t80 * t110 - t73 * t77) * t69;
t55 = t99 * t62;
t1 = [-t75 * t77 * t69, t57, t56, 0, 0, 0; (-t71 * t59 - (-t67 + (t77 * t111 + t67) * t69) * t109) * t58 (-t83 * t92 * t59 - ((t68 * t105 - t67 * t81) * t92 + t98 * t57) * t113) * t58 (t76 * t59 - (t98 * t56 - t67 * t73 + t68 * t80) * t113) * t58, 0, 0, 0; ((-t73 * t87 - t81 * t88) * t63 - (-t73 * t88 + t81 * t87) * t112) * t62 ((-t87 * t108 - t84 * t88) * t63 - (-t88 * t108 + t84 * t87) * t112) * t62 (t88 * t112 - t87 * t63) * t75 * t62, t55, t55, 0;];
Ja_rot  = t1;
