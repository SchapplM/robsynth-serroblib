% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:23
% EndTime: 2019-02-26 20:53:24
% DurationCPUTime: 0.25s
% Computational Cost: add. (976->38), mult. (2781->89), div. (55->9), fcn. (3822->17), ass. (0->59)
t103 = sin(pkin(13));
t107 = cos(pkin(13));
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t118 = t115 * t103 + t112 * t107;
t106 = sin(pkin(6));
t116 = cos(qJ(1));
t124 = t106 * t116;
t105 = sin(pkin(7));
t99 = t112 * t103 - t115 * t107;
t91 = t99 * t105;
t109 = cos(pkin(7));
t93 = t99 * t109;
t110 = cos(pkin(6));
t108 = cos(pkin(12));
t120 = t116 * t108;
t104 = sin(pkin(12));
t113 = sin(qJ(1));
t123 = t113 * t104;
t95 = -t110 * t120 + t123;
t121 = t116 * t104;
t122 = t113 * t108;
t96 = t110 * t121 + t122;
t79 = -t118 * t96 + t91 * t124 + t95 * t93;
t87 = t110 * t91 + (t104 * t118 + t108 * t93) * t106;
t73 = atan2(t79, t87);
t71 = cos(t73);
t128 = t71 * t79;
t111 = sin(qJ(5));
t114 = cos(qJ(5));
t125 = t106 * t113;
t92 = t118 * t105;
t94 = t118 * t109;
t97 = -t110 * t122 - t121;
t98 = -t110 * t123 + t120;
t83 = t92 * t125 + t97 * t94 - t98 * t99;
t89 = -t97 * t105 + t109 * t125;
t77 = t89 * t111 + t83 * t114;
t75 = 0.1e1 / t77 ^ 2;
t76 = t83 * t111 - t89 * t114;
t127 = t75 * t76;
t70 = sin(t73);
t68 = t70 * t79 + t71 * t87;
t67 = 0.1e1 / t68 ^ 2;
t81 = -t118 * t98 - t91 * t125 - t97 * t93;
t126 = t81 ^ 2 * t67;
t119 = t76 ^ 2 * t75 + 0.1e1;
t117 = t92 * t124 + t95 * t94 + t96 * t99;
t88 = -t95 * t105 + t109 * t124;
t86 = t110 * t92 + (-t104 * t99 + t108 * t94) * t106;
t85 = 0.1e1 / t87 ^ 2;
t84 = 0.1e1 / t87;
t74 = 0.1e1 / t77;
t72 = 0.1e1 / (t79 ^ 2 * t85 + 0.1e1);
t69 = 0.1e1 / t119;
t66 = 0.1e1 / t68;
t65 = 0.1e1 / (0.1e1 + t126);
t64 = (-t79 * t85 * t86 + t117 * t84) * t72;
t1 = [t81 * t84 * t72, 0, t64, 0, 0, 0; (t79 * t66 + (t70 + (t84 * t128 - t70) * t72) * t126) * t65, 0 (t83 * t66 + (t70 * t117 + t71 * t86 + (-t70 * t87 + t128) * t64) * t81 * t67) * t65, 0, 0, 0; ((t111 * t117 - t88 * t114) * t74 - (t88 * t111 + t114 * t117) * t127) * t69, 0 (t111 * t74 - t114 * t127) * t81 * t69, 0, t119 * t69, 0;];
Ja_rot  = t1;
