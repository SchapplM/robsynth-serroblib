% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:27
% EndTime: 2019-02-26 22:20:28
% DurationCPUTime: 0.32s
% Computational Cost: add. (1467->49), mult. (4186->118), div. (85->9), fcn. (5752->17), ass. (0->66)
t114 = sin(pkin(13));
t117 = cos(pkin(13));
t120 = sin(qJ(3));
t124 = cos(qJ(3));
t109 = t114 * t120 - t124 * t117;
t118 = cos(pkin(7));
t103 = t109 * t118;
t121 = sin(qJ(2));
t122 = sin(qJ(1));
t125 = cos(qJ(2));
t126 = cos(qJ(1));
t138 = cos(pkin(6));
t132 = t126 * t138;
t105 = t121 * t122 - t125 * t132;
t106 = t121 * t132 + t122 * t125;
t116 = sin(pkin(6));
t115 = sin(pkin(7));
t129 = t109 * t115;
t128 = t116 * t129;
t130 = t114 * t124 + t117 * t120;
t87 = t103 * t105 - t106 * t130 + t126 * t128;
t95 = (t103 * t125 + t121 * t130) * t116 + t138 * t129;
t81 = atan2(t87, t95);
t78 = sin(t81);
t79 = cos(t81);
t76 = t78 * t87 + t79 * t95;
t75 = 0.1e1 / t76 ^ 2;
t133 = t122 * t138;
t107 = -t126 * t121 - t125 * t133;
t108 = -t121 * t133 + t125 * t126;
t89 = -t107 * t103 - t108 * t130 - t122 * t128;
t143 = t75 * t89;
t142 = t75 * t89 ^ 2;
t141 = t79 * t87;
t136 = t116 * t122;
t100 = -t107 * t115 + t118 * t136;
t119 = sin(qJ(5));
t123 = cos(qJ(5));
t102 = t130 * t115;
t104 = t130 * t118;
t91 = t102 * t136 + t104 * t107 - t108 * t109;
t85 = t100 * t119 + t123 * t91;
t83 = 0.1e1 / t85 ^ 2;
t84 = -t100 * t123 + t119 * t91;
t140 = t83 * t84;
t93 = 0.1e1 / t95 ^ 2;
t139 = t87 * t93;
t137 = t108 * t115;
t135 = t116 * t126;
t134 = t83 * t84 ^ 2 + 0.1e1;
t131 = -t78 * t95 + t141;
t127 = t102 * t135 + t105 * t104 + t106 * t109;
t99 = -t105 * t115 + t118 * t135;
t98 = (-t103 * t121 + t125 * t130) * t116;
t97 = -t104 * t108 - t107 * t109;
t96 = t103 * t106 + t105 * t130;
t94 = t138 * t102 + (t104 * t125 - t109 * t121) * t116;
t92 = 0.1e1 / t95;
t82 = 0.1e1 / t85;
t80 = 0.1e1 / (t87 ^ 2 * t93 + 0.1e1);
t77 = 0.1e1 / t134;
t74 = 0.1e1 / t76;
t73 = 0.1e1 / (0.1e1 + t142);
t72 = (-t98 * t139 + t92 * t96) * t80;
t71 = (t127 * t92 - t94 * t139) * t80;
t1 = [t89 * t92 * t80, t72, t71, 0, 0, 0; (t87 * t74 + (t78 + (t92 * t141 - t78) * t80) * t142) * t73 ((-t103 * t108 + t107 * t130) * t74 + (t131 * t72 + t78 * t96 + t79 * t98) * t143) * t73 (t91 * t74 + (t127 * t78 + t131 * t71 + t79 * t94) * t143) * t73, 0, 0, 0; ((t119 * t127 - t123 * t99) * t82 - (t119 * t99 + t123 * t127) * t140) * t77 ((t119 * t97 - t123 * t137) * t82 - (t119 * t137 + t123 * t97) * t140) * t77 (t119 * t82 - t123 * t140) * t89 * t77, 0, t134 * t77, 0;];
Ja_rot  = t1;
