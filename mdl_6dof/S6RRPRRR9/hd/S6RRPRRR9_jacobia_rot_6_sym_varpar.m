% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR9
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
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:46
% EndTime: 2019-02-26 21:58:46
% DurationCPUTime: 0.25s
% Computational Cost: add. (1756->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
t101 = pkin(12) + qJ(4) + qJ(5);
t100 = cos(t101);
t102 = sin(pkin(6));
t109 = cos(qJ(1));
t116 = t102 * t109;
t103 = cos(pkin(6));
t105 = sin(qJ(2));
t113 = t109 * t105;
t106 = sin(qJ(1));
t108 = cos(qJ(2));
t114 = t106 * t108;
t94 = t103 * t113 + t114;
t99 = sin(t101);
t83 = t100 * t116 + t94 * t99;
t119 = t102 * t105;
t91 = -t103 * t100 + t99 * t119;
t78 = atan2(-t83, t91);
t75 = sin(t78);
t76 = cos(t78);
t73 = -t75 * t83 + t76 * t91;
t72 = 0.1e1 / t73 ^ 2;
t118 = t102 * t106;
t112 = t109 * t108;
t115 = t106 * t105;
t96 = -t103 * t115 + t112;
t87 = -t100 * t118 + t96 * t99;
t126 = t72 * t87;
t125 = t76 * t83;
t107 = cos(qJ(6));
t104 = sin(qJ(6));
t95 = t103 * t114 + t113;
t121 = t95 * t104;
t88 = t96 * t100 + t99 * t118;
t82 = t88 * t107 + t121;
t80 = 0.1e1 / t82 ^ 2;
t120 = t95 * t107;
t81 = t88 * t104 - t120;
t124 = t80 * t81;
t90 = 0.1e1 / t91 ^ 2;
t123 = t83 * t90;
t122 = t87 ^ 2 * t72;
t117 = t102 * t108;
t111 = t81 ^ 2 * t80 + 0.1e1;
t85 = t94 * t100 - t99 * t116;
t110 = -t75 * t91 - t125;
t93 = t103 * t112 - t115;
t92 = t100 * t119 + t103 * t99;
t89 = 0.1e1 / t91;
t79 = 0.1e1 / t82;
t77 = 0.1e1 / (t83 ^ 2 * t90 + 0.1e1);
t74 = 0.1e1 / t111;
t71 = 0.1e1 / t73;
t70 = 0.1e1 / (0.1e1 + t122);
t69 = (t117 * t123 - t89 * t93) * t99 * t77;
t68 = (t92 * t123 - t85 * t89) * t77;
t67 = (-t104 * t79 + t107 * t124) * t87 * t74;
t66 = (t88 * t71 - (t110 * t68 - t75 * t85 + t76 * t92) * t126) * t70;
t1 = [-t87 * t89 * t77, t69, 0, t68, t68, 0; (-t83 * t71 - (-t75 + (t89 * t125 + t75) * t77) * t122) * t70 (-t95 * t99 * t71 - ((t76 * t117 - t75 * t93) * t99 + t110 * t69) * t126) * t70, 0, t66, t66, 0; ((-t104 * t85 - t93 * t107) * t79 - (t93 * t104 - t107 * t85) * t124) * t74 ((-t100 * t121 - t96 * t107) * t79 - (-t100 * t120 + t96 * t104) * t124) * t74, 0, t67, t67, t111 * t74;];
Ja_rot  = t1;
