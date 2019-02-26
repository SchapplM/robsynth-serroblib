% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6RRRRRP9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:44:25
% EndTime: 2019-02-26 22:44:26
% DurationCPUTime: 0.20s
% Computational Cost: add. (570->38), mult. (1383->91), div. (90->9), fcn. (1970->13), ass. (0->56)
t101 = cos(qJ(1));
t94 = sin(pkin(6));
t107 = t101 * t94;
t100 = cos(qJ(2));
t98 = sin(qJ(1));
t105 = t98 * t100;
t97 = sin(qJ(2));
t106 = t101 * t97;
t95 = cos(pkin(6));
t86 = t95 * t106 + t105;
t96 = sin(qJ(3));
t99 = cos(qJ(3));
t75 = t99 * t107 + t86 * t96;
t111 = t94 * t96;
t83 = t97 * t111 - t95 * t99;
t74 = atan2(-t75, t83);
t71 = sin(t74);
t72 = cos(t74);
t65 = -t71 * t75 + t72 * t83;
t64 = 0.1e1 / t65 ^ 2;
t110 = t94 * t99;
t104 = t101 * t100;
t109 = t98 * t97;
t88 = -t95 * t109 + t104;
t79 = -t98 * t110 + t88 * t96;
t117 = t64 * t79;
t80 = t98 * t111 + t88 * t99;
t87 = t95 * t105 + t106;
t93 = qJ(4) + qJ(5);
t91 = sin(t93);
t92 = cos(t93);
t70 = t80 * t92 + t87 * t91;
t68 = 0.1e1 / t70 ^ 2;
t69 = t80 * t91 - t87 * t92;
t116 = t68 * t69;
t115 = t72 * t75;
t82 = 0.1e1 / t83 ^ 2;
t114 = t75 * t82;
t113 = t79 ^ 2 * t64;
t112 = t87 * t99;
t108 = t100 * t94;
t103 = t69 ^ 2 * t68 + 0.1e1;
t77 = -t96 * t107 + t86 * t99;
t102 = -t71 * t83 - t115;
t85 = t95 * t104 - t109;
t84 = t97 * t110 + t95 * t96;
t81 = 0.1e1 / t83;
t73 = 0.1e1 / (t75 ^ 2 * t82 + 0.1e1);
t67 = 0.1e1 / t70;
t66 = 0.1e1 / t103;
t63 = 0.1e1 / t65;
t62 = 0.1e1 / (0.1e1 + t113);
t61 = (t108 * t114 - t81 * t85) * t96 * t73;
t60 = (t84 * t114 - t77 * t81) * t73;
t59 = t103 * t66;
t1 = [-t79 * t81 * t73, t61, t60, 0, 0, 0; (-t75 * t63 - (-t71 + (t81 * t115 + t71) * t73) * t113) * t62 (-t87 * t96 * t63 - ((t72 * t108 - t71 * t85) * t96 + t102 * t61) * t117) * t62 (t80 * t63 - (t102 * t60 - t71 * t77 + t72 * t84) * t117) * t62, 0, 0, 0; ((-t77 * t91 - t85 * t92) * t67 - (-t77 * t92 + t85 * t91) * t116) * t66 ((-t91 * t112 - t88 * t92) * t67 - (-t92 * t112 + t88 * t91) * t116) * t66 (t92 * t116 - t91 * t67) * t79 * t66, t59, t59, 0;];
Ja_rot  = t1;
