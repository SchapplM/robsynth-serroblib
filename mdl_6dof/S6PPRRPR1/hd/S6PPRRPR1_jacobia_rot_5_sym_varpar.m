% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRPR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:29
% EndTime: 2019-02-26 19:40:29
% DurationCPUTime: 0.23s
% Computational Cost: add. (961->43), mult. (2761->105), div. (60->9), fcn. (3790->17), ass. (0->64)
t89 = sin(pkin(11));
t96 = cos(pkin(6));
t114 = t89 * t96;
t88 = sin(pkin(12));
t93 = cos(pkin(12));
t94 = cos(pkin(11));
t105 = t93 * t114 + t94 * t88;
t90 = sin(pkin(7));
t91 = sin(pkin(6));
t112 = t91 * t90;
t95 = cos(pkin(7));
t119 = t105 * t95 - t89 * t112;
t109 = t94 * t96;
t106 = t93 * t109 - t89 * t88;
t111 = t91 * t95;
t103 = -t106 * t90 - t94 * t111;
t100 = cos(qJ(3));
t102 = t106 * t95 - t94 * t112;
t84 = t88 * t109 + t89 * t93;
t98 = sin(qJ(3));
t71 = t84 * t100 + t102 * t98;
t97 = sin(qJ(4));
t99 = cos(qJ(4));
t65 = -t103 * t99 + t71 * t97;
t110 = t93 * t95;
t113 = t90 * t96;
t81 = t98 * t113 + (t100 * t88 + t98 * t110) * t91;
t83 = -t93 * t112 + t96 * t95;
t76 = t81 * t97 - t83 * t99;
t64 = atan2(-t65, t76);
t61 = sin(t64);
t62 = cos(t64);
t55 = -t61 * t65 + t62 * t76;
t54 = 0.1e1 / t55 ^ 2;
t101 = t105 * t90 + t89 * t111;
t85 = -t88 * t114 + t94 * t93;
t73 = t85 * t100 - t119 * t98;
t68 = -t101 * t99 + t73 * t97;
t118 = t54 * t68;
t69 = t101 * t97 + t73 * t99;
t72 = t119 * t100 + t85 * t98;
t87 = sin(pkin(13));
t92 = cos(pkin(13));
t60 = t69 * t92 + t72 * t87;
t58 = 0.1e1 / t60 ^ 2;
t59 = t69 * t87 - t72 * t92;
t117 = t58 * t59;
t75 = 0.1e1 / t76 ^ 2;
t116 = t65 * t75;
t115 = t72 * t99;
t107 = -t61 * t76 - t62 * t65;
t80 = -t91 * t88 * t98 + (t91 * t110 + t113) * t100;
t77 = t81 * t99 + t83 * t97;
t74 = 0.1e1 / t76;
t70 = t102 * t100 - t84 * t98;
t67 = t103 * t97 + t71 * t99;
t63 = 0.1e1 / (t65 ^ 2 * t75 + 0.1e1);
t57 = 0.1e1 / t60;
t56 = 0.1e1 / (t59 ^ 2 * t58 + 0.1e1);
t53 = 0.1e1 / t55;
t52 = 0.1e1 / (t68 ^ 2 * t54 + 0.1e1);
t51 = (t80 * t116 - t70 * t74) * t97 * t63;
t50 = (t77 * t116 - t67 * t74) * t63;
t1 = [0, 0, t51, t50, 0, 0; 0, 0 (-t72 * t97 * t53 - ((-t61 * t70 + t62 * t80) * t97 + t107 * t51) * t118) * t52 (t69 * t53 - (t107 * t50 - t61 * t67 + t62 * t77) * t118) * t52, 0, 0; 0, 0 ((-t87 * t115 - t73 * t92) * t57 - (-t92 * t115 + t73 * t87) * t117) * t56 (t92 * t117 - t57 * t87) * t68 * t56, 0, 0;];
Ja_rot  = t1;
