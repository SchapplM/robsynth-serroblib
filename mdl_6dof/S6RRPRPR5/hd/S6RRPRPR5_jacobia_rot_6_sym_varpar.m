% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:08
% EndTime: 2019-02-26 21:40:08
% DurationCPUTime: 0.27s
% Computational Cost: add. (933->40), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->58)
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t104 = sin(pkin(6));
t112 = cos(qJ(1));
t118 = t104 * t112;
t109 = sin(qJ(1));
t106 = cos(pkin(6));
t103 = sin(pkin(11));
t105 = cos(pkin(11));
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t114 = t111 * t103 + t108 * t105;
t94 = t114 * t106;
t95 = t108 * t103 - t111 * t105;
t83 = -t109 * t95 + t112 * t94;
t76 = t83 * t107 + t110 * t118;
t93 = t114 * t104;
t89 = -t106 * t110 + t93 * t107;
t75 = atan2(-t76, t89);
t72 = sin(t75);
t73 = cos(t75);
t66 = -t72 * t76 + t73 * t89;
t65 = 0.1e1 / t66 ^ 2;
t115 = -t109 * t94 - t112 * t95;
t119 = t104 * t109;
t80 = t107 * t115 - t110 * t119;
t126 = t65 * t80;
t102 = pkin(12) + qJ(6);
t101 = cos(t102);
t100 = sin(t102);
t113 = t95 * t106;
t85 = t109 * t113 - t112 * t114;
t121 = t85 * t100;
t81 = t107 * t119 + t110 * t115;
t71 = t81 * t101 - t121;
t69 = 0.1e1 / t71 ^ 2;
t120 = t85 * t101;
t70 = t81 * t100 + t120;
t125 = t69 * t70;
t124 = t73 * t76;
t88 = 0.1e1 / t89 ^ 2;
t123 = t76 * t88;
t122 = t80 ^ 2 * t65;
t117 = t70 ^ 2 * t69 + 0.1e1;
t78 = -t107 * t118 + t83 * t110;
t116 = -t72 * t89 - t124;
t92 = t95 * t104;
t90 = t106 * t107 + t93 * t110;
t87 = 0.1e1 / t89;
t82 = -t109 * t114 - t112 * t113;
t74 = 0.1e1 / (t76 ^ 2 * t88 + 0.1e1);
t68 = 0.1e1 / t71;
t67 = 0.1e1 / t117;
t64 = 0.1e1 / t66;
t63 = 0.1e1 / (0.1e1 + t122);
t62 = (-t123 * t92 - t82 * t87) * t74 * t107;
t61 = (t123 * t90 - t78 * t87) * t74;
t1 = [-t80 * t87 * t74, t62, 0, t61, 0, 0; (-t76 * t64 - (-t72 + (t124 * t87 + t72) * t74) * t122) * t63 (t85 * t107 * t64 - (t116 * t62 + (-t72 * t82 - t73 * t92) * t107) * t126) * t63, 0 (t81 * t64 - (t116 * t61 - t72 * t78 + t73 * t90) * t126) * t63, 0, 0; ((-t100 * t78 - t82 * t101) * t68 - (t82 * t100 - t101 * t78) * t125) * t67 ((-t101 * t115 + t110 * t121) * t68 - (t100 * t115 + t110 * t120) * t125) * t67, 0 (-t100 * t68 + t101 * t125) * t80 * t67, 0, t117 * t67;];
Ja_rot  = t1;
