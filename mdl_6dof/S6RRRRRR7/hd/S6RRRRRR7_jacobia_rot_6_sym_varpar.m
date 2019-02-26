% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:41
% EndTime: 2019-02-26 22:50:42
% DurationCPUTime: 0.24s
% Computational Cost: add. (709->38), mult. (1466->91), div. (95->9), fcn. (2082->13), ass. (0->56)
t100 = cos(qJ(3));
t102 = cos(qJ(1));
t95 = sin(pkin(6));
t108 = t102 * t95;
t101 = cos(qJ(2));
t99 = sin(qJ(1));
t106 = t99 * t101;
t98 = sin(qJ(2));
t107 = t102 * t98;
t96 = cos(pkin(6));
t87 = t96 * t107 + t106;
t97 = sin(qJ(3));
t76 = t100 * t108 + t87 * t97;
t113 = t95 * t97;
t84 = -t96 * t100 + t98 * t113;
t75 = atan2(-t76, t84);
t72 = sin(t75);
t73 = cos(t75);
t66 = -t72 * t76 + t73 * t84;
t65 = 0.1e1 / t66 ^ 2;
t110 = t100 * t95;
t105 = t102 * t101;
t112 = t99 * t98;
t89 = -t96 * t112 + t105;
t80 = -t99 * t110 + t89 * t97;
t118 = t65 * t80;
t81 = t89 * t100 + t99 * t113;
t88 = t96 * t106 + t107;
t94 = qJ(4) + qJ(5) + qJ(6);
t92 = sin(t94);
t93 = cos(t94);
t71 = t81 * t93 + t88 * t92;
t69 = 0.1e1 / t71 ^ 2;
t70 = t81 * t92 - t88 * t93;
t117 = t69 * t70;
t116 = t73 * t76;
t83 = 0.1e1 / t84 ^ 2;
t115 = t76 * t83;
t114 = t80 ^ 2 * t65;
t111 = t100 * t88;
t109 = t101 * t95;
t104 = t70 ^ 2 * t69 + 0.1e1;
t78 = t87 * t100 - t97 * t108;
t103 = -t72 * t84 - t116;
t86 = t96 * t105 - t112;
t85 = t98 * t110 + t96 * t97;
t82 = 0.1e1 / t84;
t74 = 0.1e1 / (t76 ^ 2 * t83 + 0.1e1);
t68 = 0.1e1 / t71;
t67 = 0.1e1 / t104;
t64 = 0.1e1 / t66;
t63 = 0.1e1 / (0.1e1 + t114);
t62 = (t109 * t115 - t82 * t86) * t97 * t74;
t61 = (t85 * t115 - t78 * t82) * t74;
t60 = t104 * t67;
t1 = [-t80 * t82 * t74, t62, t61, 0, 0, 0; (-t76 * t64 - (-t72 + (t82 * t116 + t72) * t74) * t114) * t63 (-t88 * t97 * t64 - ((t73 * t109 - t72 * t86) * t97 + t103 * t62) * t118) * t63 (t81 * t64 - (t103 * t61 - t72 * t78 + t73 * t85) * t118) * t63, 0, 0, 0; ((-t78 * t92 - t86 * t93) * t68 - (-t78 * t93 + t86 * t92) * t117) * t67 ((-t92 * t111 - t89 * t93) * t68 - (-t93 * t111 + t89 * t92) * t117) * t67 (t93 * t117 - t92 * t68) * t80 * t67, t60, t60, t60;];
Ja_rot  = t1;
