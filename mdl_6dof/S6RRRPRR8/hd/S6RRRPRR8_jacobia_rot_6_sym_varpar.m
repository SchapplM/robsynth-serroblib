% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6RRRPRR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:57
% EndTime: 2019-02-26 22:19:57
% DurationCPUTime: 0.23s
% Computational Cost: add. (1003->39), mult. (1383->90), div. (90->9), fcn. (1970->13), ass. (0->58)
t104 = cos(qJ(1));
t99 = sin(pkin(6));
t111 = t104 * t99;
t100 = cos(pkin(6));
t101 = sin(qJ(2));
t108 = t104 * t101;
t102 = sin(qJ(1));
t103 = cos(qJ(2));
t109 = t102 * t103;
t88 = t100 * t108 + t109;
t97 = qJ(3) + pkin(12);
t93 = sin(t97);
t94 = cos(t97);
t77 = t94 * t111 + t88 * t93;
t114 = t101 * t99;
t85 = -t100 * t94 + t93 * t114;
t76 = atan2(-t77, t85);
t73 = sin(t76);
t74 = cos(t76);
t67 = -t73 * t77 + t74 * t85;
t66 = 0.1e1 / t67 ^ 2;
t113 = t102 * t99;
t107 = t104 * t103;
t110 = t102 * t101;
t90 = -t100 * t110 + t107;
t81 = -t94 * t113 + t90 * t93;
t121 = t66 * t81;
t89 = t100 * t109 + t108;
t98 = qJ(5) + qJ(6);
t95 = sin(t98);
t116 = t89 * t95;
t82 = t93 * t113 + t90 * t94;
t96 = cos(t98);
t72 = t82 * t96 + t116;
t70 = 0.1e1 / t72 ^ 2;
t115 = t89 * t96;
t71 = t82 * t95 - t115;
t120 = t70 * t71;
t119 = t74 * t77;
t84 = 0.1e1 / t85 ^ 2;
t118 = t77 * t84;
t117 = t81 ^ 2 * t66;
t112 = t103 * t99;
t106 = t71 ^ 2 * t70 + 0.1e1;
t79 = -t93 * t111 + t88 * t94;
t105 = -t73 * t85 - t119;
t87 = t100 * t107 - t110;
t86 = t100 * t93 + t94 * t114;
t83 = 0.1e1 / t85;
t75 = 0.1e1 / (t77 ^ 2 * t84 + 0.1e1);
t69 = 0.1e1 / t72;
t68 = 0.1e1 / t106;
t65 = 0.1e1 / t67;
t64 = 0.1e1 / (0.1e1 + t117);
t63 = (t112 * t118 - t83 * t87) * t93 * t75;
t62 = (t86 * t118 - t79 * t83) * t75;
t61 = t106 * t68;
t1 = [-t81 * t83 * t75, t63, t62, 0, 0, 0; (-t77 * t65 - (-t73 + (t83 * t119 + t73) * t75) * t117) * t64 (-t89 * t93 * t65 - ((t74 * t112 - t73 * t87) * t93 + t105 * t63) * t121) * t64 (t82 * t65 - (t105 * t62 - t73 * t79 + t74 * t86) * t121) * t64, 0, 0, 0; ((-t79 * t95 - t87 * t96) * t69 - (-t79 * t96 + t87 * t95) * t120) * t68 ((-t94 * t116 - t90 * t96) * t69 - (-t94 * t115 + t90 * t95) * t120) * t68 (t96 * t120 - t95 * t69) * t81 * t68, 0, t61, t61;];
Ja_rot  = t1;
