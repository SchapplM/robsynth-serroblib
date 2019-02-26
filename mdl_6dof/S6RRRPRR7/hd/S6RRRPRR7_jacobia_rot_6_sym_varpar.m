% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR7
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
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:18
% EndTime: 2019-02-26 22:19:19
% DurationCPUTime: 0.25s
% Computational Cost: add. (1756->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
t102 = qJ(3) + pkin(12) + qJ(5);
t100 = sin(t102);
t101 = cos(t102);
t103 = sin(pkin(6));
t110 = cos(qJ(1));
t117 = t103 * t110;
t104 = cos(pkin(6));
t106 = sin(qJ(2));
t114 = t110 * t106;
t107 = sin(qJ(1));
t109 = cos(qJ(2));
t115 = t107 * t109;
t95 = t104 * t114 + t115;
t84 = t95 * t100 + t101 * t117;
t120 = t103 * t106;
t92 = t100 * t120 - t104 * t101;
t79 = atan2(-t84, t92);
t76 = sin(t79);
t77 = cos(t79);
t74 = -t76 * t84 + t77 * t92;
t73 = 0.1e1 / t74 ^ 2;
t119 = t103 * t107;
t113 = t110 * t109;
t116 = t107 * t106;
t97 = -t104 * t116 + t113;
t88 = t97 * t100 - t101 * t119;
t127 = t73 * t88;
t126 = t77 * t84;
t108 = cos(qJ(6));
t105 = sin(qJ(6));
t96 = t104 * t115 + t114;
t122 = t96 * t105;
t89 = t100 * t119 + t97 * t101;
t83 = t89 * t108 + t122;
t81 = 0.1e1 / t83 ^ 2;
t121 = t96 * t108;
t82 = t89 * t105 - t121;
t125 = t81 * t82;
t91 = 0.1e1 / t92 ^ 2;
t124 = t84 * t91;
t123 = t88 ^ 2 * t73;
t118 = t103 * t109;
t112 = t82 ^ 2 * t81 + 0.1e1;
t86 = -t100 * t117 + t95 * t101;
t111 = -t76 * t92 - t126;
t94 = t104 * t113 - t116;
t93 = t104 * t100 + t101 * t120;
t90 = 0.1e1 / t92;
t80 = 0.1e1 / t83;
t78 = 0.1e1 / (t84 ^ 2 * t91 + 0.1e1);
t75 = 0.1e1 / t112;
t72 = 0.1e1 / t74;
t71 = 0.1e1 / (0.1e1 + t123);
t70 = (t118 * t124 - t90 * t94) * t78 * t100;
t69 = (t93 * t124 - t86 * t90) * t78;
t68 = (-t105 * t80 + t108 * t125) * t88 * t75;
t67 = (t89 * t72 - (t111 * t69 - t76 * t86 + t77 * t93) * t127) * t71;
t1 = [-t88 * t90 * t78, t70, t69, 0, t69, 0; (-t84 * t72 - (-t76 + (t90 * t126 + t76) * t78) * t123) * t71 (-t96 * t100 * t72 - (t111 * t70 + (t77 * t118 - t76 * t94) * t100) * t127) * t71, t67, 0, t67, 0; ((-t105 * t86 - t94 * t108) * t80 - (t94 * t105 - t108 * t86) * t125) * t75 ((-t101 * t122 - t97 * t108) * t80 - (-t101 * t121 + t97 * t105) * t125) * t75, t68, 0, t68, t112 * t75;];
Ja_rot  = t1;
