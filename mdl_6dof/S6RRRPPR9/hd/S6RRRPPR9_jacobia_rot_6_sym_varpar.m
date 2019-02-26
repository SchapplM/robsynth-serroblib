% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:09
% EndTime: 2019-02-26 22:08:10
% DurationCPUTime: 0.23s
% Computational Cost: add. (623->43), mult. (1730->106), div. (85->9), fcn. (2452->15), ass. (0->63)
t129 = sin(qJ(1));
t105 = sin(qJ(3));
t108 = cos(qJ(3));
t110 = cos(qJ(1));
t121 = sin(pkin(6));
t113 = t110 * t121;
t103 = cos(pkin(6));
t109 = cos(qJ(2));
t116 = t129 * t109;
t106 = sin(qJ(2));
t120 = t110 * t106;
t96 = t103 * t120 + t116;
t85 = t96 * t105 + t108 * t113;
t115 = t106 * t121;
t93 = t103 * t108 - t105 * t115;
t82 = atan2(t85, t93);
t79 = sin(t82);
t80 = cos(t82);
t70 = t79 * t85 + t80 * t93;
t69 = 0.1e1 / t70 ^ 2;
t112 = t121 * t129;
t117 = t129 * t106;
t119 = t110 * t109;
t98 = -t103 * t117 + t119;
t88 = t98 * t105 - t108 * t112;
t128 = t69 * t88;
t104 = sin(qJ(6));
t107 = cos(qJ(6));
t101 = sin(pkin(11));
t102 = cos(pkin(11));
t97 = t103 * t116 + t120;
t122 = t97 * t102;
t90 = t105 * t112 + t98 * t108;
t77 = t90 * t101 - t122;
t123 = t97 * t101;
t78 = t90 * t102 + t123;
t74 = t77 * t104 + t78 * t107;
t72 = 0.1e1 / t74 ^ 2;
t73 = t78 * t104 - t77 * t107;
t127 = t72 * t73;
t126 = t80 * t85;
t92 = 0.1e1 / t93 ^ 2;
t125 = t85 * t92;
t124 = t88 ^ 2 * t69;
t118 = t73 ^ 2 * t72 + 0.1e1;
t114 = t109 * t121;
t86 = -t105 * t113 + t96 * t108;
t111 = -t79 * t93 + t126;
t95 = t103 * t119 - t117;
t94 = -t103 * t105 - t108 * t115;
t91 = 0.1e1 / t93;
t84 = t98 * t101 - t108 * t122;
t83 = -t98 * t102 - t108 * t123;
t81 = 0.1e1 / (t85 ^ 2 * t92 + 0.1e1);
t76 = t95 * t101 - t102 * t86;
t75 = -t101 * t86 - t95 * t102;
t71 = 0.1e1 / t74;
t68 = 0.1e1 / t70;
t67 = 0.1e1 / (0.1e1 + t124);
t66 = (t114 * t125 + t91 * t95) * t81 * t105;
t65 = (-t94 * t125 + t86 * t91) * t81;
t64 = 0.1e1 / t118;
t1 = [t88 * t91 * t81, t66, t65, 0, 0, 0; (t85 * t68 + (t79 + (t91 * t126 - t79) * t81) * t124) * t67 (t97 * t105 * t68 + (t111 * t66 + (-t80 * t114 + t79 * t95) * t105) * t128) * t67 (-t90 * t68 + (t111 * t65 + t79 * t86 + t80 * t94) * t128) * t67, 0, 0, 0; ((t76 * t104 - t75 * t107) * t71 - (t75 * t104 + t76 * t107) * t127) * t64 ((t84 * t104 - t83 * t107) * t71 - (t83 * t104 + t84 * t107) * t127) * t64 ((t101 * t107 - t102 * t104) * t71 - (-t101 * t104 - t102 * t107) * t127) * t64 * t88, 0, 0, t118 * t64;];
Ja_rot  = t1;
