% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:47
% EndTime: 2019-02-26 20:01:48
% DurationCPUTime: 0.26s
% Computational Cost: add. (1149->41), mult. (1993->102), div. (87->9), fcn. (2807->13), ass. (0->60)
t101 = sin(qJ(5));
t102 = sin(qJ(2));
t100 = cos(pkin(6));
t104 = cos(qJ(2));
t108 = t100 * t104;
t97 = sin(pkin(10));
t99 = cos(pkin(10));
t105 = -t97 * t102 + t99 * t108;
t119 = t105 * t101;
t103 = cos(qJ(5));
t98 = sin(pkin(6));
t113 = t98 * t99;
t109 = t100 * t102;
t91 = t97 * t104 + t99 * t109;
t96 = qJ(3) + pkin(11);
t94 = sin(t96);
t95 = cos(t96);
t80 = -t94 * t113 + t91 * t95;
t72 = t80 * t101 + t105 * t103;
t112 = t102 * t98;
t89 = t100 * t94 + t95 * t112;
t85 = t98 * t104 * t103 + t89 * t101;
t69 = atan2(-t72, t85);
t66 = sin(t69);
t67 = cos(t69);
t64 = -t66 * t72 + t67 * t85;
t63 = 0.1e1 / t64 ^ 2;
t92 = t99 * t102 + t97 * t108;
t110 = t92 * t103;
t114 = t97 * t98;
t93 = t99 * t104 - t97 * t109;
t82 = t94 * t114 + t93 * t95;
t75 = t82 * t101 - t110;
t118 = t63 * t75;
t111 = t92 * t101;
t76 = t82 * t103 + t111;
t71 = 0.1e1 / t76 ^ 2;
t81 = t95 * t114 - t93 * t94;
t117 = t71 * t81;
t84 = 0.1e1 / t85 ^ 2;
t116 = t72 * t84;
t115 = t81 ^ 2 * t71;
t107 = t101 * t104;
t106 = -t66 * t85 - t67 * t72;
t88 = t100 * t95 - t94 * t112;
t87 = (-t102 * t103 + t95 * t107) * t98;
t86 = t89 * t103 - t98 * t107;
t83 = 0.1e1 / t85;
t79 = -t95 * t113 - t91 * t94;
t77 = -t91 * t103 + t95 * t119;
t74 = t80 * t103 - t119;
t70 = 0.1e1 / t76;
t68 = 0.1e1 / (t72 ^ 2 * t84 + 0.1e1);
t65 = 0.1e1 / (0.1e1 + t115);
t62 = 0.1e1 / t64;
t61 = 0.1e1 / (t75 ^ 2 * t63 + 0.1e1);
t60 = (t88 * t116 - t79 * t83) * t68 * t101;
t59 = (t87 * t116 - t77 * t83) * t68;
t58 = (t86 * t116 - t74 * t83) * t68;
t1 = [0, t59, t60, 0, t58, 0; 0 ((-t93 * t103 - t95 * t111) * t62 - (t106 * t59 - t66 * t77 + t67 * t87) * t118) * t61 (t81 * t101 * t62 - (t106 * t60 + (-t66 * t79 + t67 * t88) * t101) * t118) * t61, 0 (t76 * t62 - (t106 * t58 - t66 * t74 + t67 * t86) * t118) * t61, 0; 0 (t92 * t94 * t70 - (t93 * t101 - t95 * t110) * t117) * t65 (-t103 * t115 - t70 * t82) * t65, 0, t75 * t65 * t117, 0;];
Ja_rot  = t1;
