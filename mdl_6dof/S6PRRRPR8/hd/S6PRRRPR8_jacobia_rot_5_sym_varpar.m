% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR8_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:42
% EndTime: 2019-02-26 20:14:42
% DurationCPUTime: 0.26s
% Computational Cost: add. (1245->53), mult. (3597->128), div. (87->9), fcn. (4945->15), ass. (0->67)
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t109 = sin(pkin(7));
t111 = cos(pkin(12));
t108 = sin(pkin(12));
t116 = sin(qJ(2));
t113 = cos(pkin(6));
t119 = cos(qJ(2));
t127 = t113 * t119;
t122 = -t108 * t116 + t111 * t127;
t110 = sin(pkin(6));
t112 = cos(pkin(7));
t131 = t110 * t112;
t121 = -t122 * t109 - t111 * t131;
t128 = t113 * t116;
t103 = t108 * t119 + t111 * t128;
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t132 = t110 * t109;
t120 = -t111 * t132 + t122 * t112;
t90 = t103 * t118 + t120 * t115;
t79 = t90 * t114 - t121 * t117;
t102 = t113 * t112 - t119 * t132;
t130 = t112 * t115;
t134 = t109 * t113;
t99 = t115 * t134 + (t116 * t118 + t119 * t130) * t110;
t93 = -t102 * t117 + t99 * t114;
t78 = atan2(-t79, t93);
t74 = sin(t78);
t75 = cos(t78);
t73 = -t74 * t79 + t75 * t93;
t72 = 0.1e1 / t73 ^ 2;
t104 = -t108 * t127 - t111 * t116;
t100 = -t104 * t109 + t108 * t131;
t105 = -t108 * t128 + t111 * t119;
t124 = t108 * t132;
t92 = t105 * t118 + (t104 * t112 + t124) * t115;
t82 = -t100 * t117 + t92 * t114;
t137 = t72 * t82;
t88 = 0.1e1 / t93 ^ 2;
t136 = t79 * t88;
t83 = t100 * t114 + t92 * t117;
t129 = t112 * t118;
t91 = -t104 * t129 + t105 * t115 - t118 * t124;
t86 = 0.1e1 / t91 ^ 2;
t135 = t83 * t86;
t133 = t109 * t117;
t126 = t115 * t116;
t125 = t118 * t119;
t123 = -t74 * t93 - t75 * t79;
t98 = t118 * t134 + (t112 * t125 - t126) * t110;
t96 = ((-t112 * t126 + t125) * t114 - t116 * t133) * t110;
t95 = t104 * t118 - t105 * t130;
t94 = t102 * t114 + t99 * t117;
t89 = -t103 * t115 + t120 * t118;
t87 = 0.1e1 / t93;
t85 = 0.1e1 / t91;
t84 = (-t103 * t130 + t122 * t118) * t114 - t103 * t133;
t81 = t121 * t114 + t90 * t117;
t77 = 0.1e1 / (t79 ^ 2 * t88 + 0.1e1);
t76 = 0.1e1 / (t83 ^ 2 * t86 + 0.1e1);
t71 = 0.1e1 / t73;
t70 = 0.1e1 / (t82 ^ 2 * t72 + 0.1e1);
t69 = (t98 * t136 - t87 * t89) * t77 * t114;
t68 = (t96 * t136 - t84 * t87) * t77;
t67 = (t94 * t136 - t81 * t87) * t77;
t1 = [0, t68, t69, t67, 0, 0; 0 ((-t105 * t133 + t95 * t114) * t71 - (t123 * t68 - t74 * t84 + t75 * t96) * t137) * t70 (-t91 * t114 * t71 - (t123 * t69 + (-t74 * t89 + t75 * t98) * t114) * t137) * t70 (t83 * t71 - (t123 * t67 - t74 * t81 + t75 * t94) * t137) * t70, 0, 0; 0 ((t105 * t109 * t114 + t95 * t117) * t85 - (t104 * t115 + t105 * t129) * t135) * t76 (-t117 * t85 * t91 - t92 * t135) * t76, -t82 * t85 * t76, 0, 0;];
Ja_rot  = t1;
