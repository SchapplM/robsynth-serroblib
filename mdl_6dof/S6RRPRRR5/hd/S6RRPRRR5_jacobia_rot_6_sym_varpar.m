% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:19
% EndTime: 2019-02-26 21:56:19
% DurationCPUTime: 0.26s
% Computational Cost: add. (1006->40), mult. (2516->96), div. (90->9), fcn. (3539->15), ass. (0->59)
t115 = sin(pkin(6));
t114 = sin(pkin(12));
t116 = cos(pkin(12));
t119 = sin(qJ(2));
t122 = cos(qJ(2));
t125 = t122 * t114 + t119 * t116;
t104 = t125 * t115;
t117 = cos(pkin(6));
t118 = sin(qJ(4));
t121 = cos(qJ(4));
t100 = t104 * t118 - t117 * t121;
t123 = cos(qJ(1));
t129 = t115 * t123;
t105 = t125 * t117;
t106 = t119 * t114 - t122 * t116;
t120 = sin(qJ(1));
t94 = t123 * t105 - t120 * t106;
t87 = t94 * t118 + t121 * t129;
t86 = atan2(-t87, t100);
t83 = sin(t86);
t84 = cos(t86);
t77 = t84 * t100 - t83 * t87;
t76 = 0.1e1 / t77 ^ 2;
t126 = -t120 * t105 - t123 * t106;
t130 = t115 * t120;
t91 = t118 * t126 - t121 * t130;
t137 = t76 * t91;
t113 = qJ(5) + qJ(6);
t112 = cos(t113);
t111 = sin(t113);
t124 = t106 * t117;
t96 = t120 * t124 - t123 * t125;
t132 = t96 * t111;
t92 = t118 * t130 + t121 * t126;
t82 = t92 * t112 - t132;
t80 = 0.1e1 / t82 ^ 2;
t131 = t96 * t112;
t81 = t92 * t111 + t131;
t136 = t80 * t81;
t135 = t84 * t87;
t99 = 0.1e1 / t100 ^ 2;
t134 = t87 * t99;
t133 = t91 ^ 2 * t76;
t128 = t81 ^ 2 * t80 + 0.1e1;
t89 = -t118 * t129 + t94 * t121;
t127 = -t100 * t83 - t135;
t103 = t106 * t115;
t101 = t104 * t121 + t117 * t118;
t98 = 0.1e1 / t100;
t93 = -t120 * t125 - t123 * t124;
t85 = 0.1e1 / (t87 ^ 2 * t99 + 0.1e1);
t79 = 0.1e1 / t82;
t78 = 0.1e1 / t128;
t75 = 0.1e1 / t77;
t74 = 0.1e1 / (0.1e1 + t133);
t73 = (-t103 * t134 - t93 * t98) * t85 * t118;
t72 = (t101 * t134 - t89 * t98) * t85;
t71 = t128 * t78;
t1 = [-t91 * t98 * t85, t73, 0, t72, 0, 0; (-t87 * t75 - (-t83 + (t98 * t135 + t83) * t85) * t133) * t74 (t96 * t118 * t75 - (t127 * t73 + (-t103 * t84 - t83 * t93) * t118) * t137) * t74, 0 (t92 * t75 - (t84 * t101 + t127 * t72 - t83 * t89) * t137) * t74, 0, 0; ((-t111 * t89 - t93 * t112) * t79 - (t93 * t111 - t112 * t89) * t136) * t78 ((-t112 * t126 + t121 * t132) * t79 - (t111 * t126 + t121 * t131) * t136) * t78, 0 (-t111 * t79 + t112 * t136) * t91 * t78, t71, t71;];
Ja_rot  = t1;
