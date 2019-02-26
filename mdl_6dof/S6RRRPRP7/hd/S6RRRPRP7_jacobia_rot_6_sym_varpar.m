% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:38
% EndTime: 2019-02-26 22:12:39
% DurationCPUTime: 0.33s
% Computational Cost: add. (1452->45), mult. (2519->109), div. (107->9), fcn. (3557->13), ass. (0->64)
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t111 = cos(pkin(6));
t116 = cos(qJ(2));
t117 = cos(qJ(1));
t122 = t117 * t116;
t113 = sin(qJ(2));
t114 = sin(qJ(1));
t125 = t114 * t113;
t121 = -t111 * t122 + t125;
t123 = t117 * t113;
t124 = t114 * t116;
t103 = t111 * t123 + t124;
t109 = qJ(3) + pkin(11);
t107 = sin(t109);
t108 = cos(t109);
t110 = sin(pkin(6));
t127 = t110 * t117;
t95 = -t103 * t108 + t107 * t127;
t140 = t95 * t112 + t121 * t115;
t119 = t121 * t112;
t139 = t95 * t115 - t119;
t129 = t110 * t113;
t100 = t111 * t107 + t108 * t129;
t91 = t110 * t116 * t115 + t100 * t112;
t79 = atan2(t140, t91);
t75 = sin(t79);
t76 = cos(t79);
t74 = t140 * t75 + t76 * t91;
t73 = 0.1e1 / t74 ^ 2;
t104 = t111 * t124 + t123;
t130 = t104 * t115;
t105 = -t111 * t125 + t122;
t128 = t110 * t114;
t97 = t105 * t108 + t107 * t128;
t85 = t97 * t112 - t130;
t137 = t73 * t85;
t136 = t76 * t140;
t131 = t104 * t112;
t86 = t97 * t115 + t131;
t81 = 0.1e1 / t86 ^ 2;
t96 = -t105 * t107 + t108 * t128;
t135 = t81 * t96;
t89 = 0.1e1 / t91 ^ 2;
t134 = t140 * t89;
t133 = t85 ^ 2 * t73;
t132 = t96 ^ 2 * t81;
t126 = t112 * t116;
t120 = -t75 * t91 + t136;
t118 = t103 * t107 + t108 * t127;
t99 = -t107 * t129 + t111 * t108;
t98 = (t108 * t126 - t113 * t115) * t110;
t92 = t100 * t115 - t110 * t126;
t88 = 0.1e1 / t91;
t87 = -t103 * t115 - t108 * t119;
t80 = 0.1e1 / t86;
t78 = 0.1e1 / (0.1e1 + t132);
t77 = 0.1e1 / (t140 ^ 2 * t89 + 0.1e1);
t72 = 0.1e1 / t74;
t71 = 0.1e1 / (0.1e1 + t133);
t70 = (t118 * t88 - t99 * t134) * t77 * t112;
t69 = (-t98 * t134 - t87 * t88) * t77;
t68 = (-t92 * t134 + t139 * t88) * t77;
t1 = [-t85 * t88 * t77, t69, t70, 0, t68, 0; (t140 * t72 - (-t75 + (-t88 * t136 + t75) * t77) * t133) * t71 ((-t105 * t115 - t108 * t131) * t72 - (t120 * t69 - t75 * t87 + t76 * t98) * t137) * t71 (t96 * t112 * t72 - (t120 * t70 + (t118 * t75 + t76 * t99) * t112) * t137) * t71, 0 (t86 * t72 - (t120 * t68 + t139 * t75 + t76 * t92) * t137) * t71, 0; (t118 * t80 - t139 * t135) * t78 (t104 * t107 * t80 - (t105 * t112 - t108 * t130) * t135) * t78 (-t115 * t132 - t80 * t97) * t78, 0, t85 * t78 * t135, 0;];
Ja_rot  = t1;
