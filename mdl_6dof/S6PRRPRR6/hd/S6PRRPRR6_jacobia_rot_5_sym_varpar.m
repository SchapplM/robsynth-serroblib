% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:05
% EndTime: 2019-02-26 20:07:05
% DurationCPUTime: 0.21s
% Computational Cost: add. (655->41), mult. (1825->98), div. (65->9), fcn. (2498->15), ass. (0->61)
t100 = cos(pkin(12));
t103 = sin(qJ(3));
t105 = cos(qJ(3));
t101 = cos(pkin(7));
t104 = sin(qJ(2));
t102 = cos(pkin(6));
t106 = cos(qJ(2));
t114 = t102 * t106;
t97 = sin(pkin(12));
t108 = t100 * t114 - t97 * t104;
t107 = t108 * t101;
t98 = sin(pkin(7));
t99 = sin(pkin(6));
t119 = t98 * t99;
t112 = t105 * t119;
t115 = t102 * t104;
t89 = t100 * t115 + t97 * t106;
t74 = t100 * t112 + t89 * t103 - t105 * t107;
t118 = t101 * t99;
t109 = t102 * t98 + t106 * t118;
t83 = t99 * t104 * t103 - t109 * t105;
t73 = atan2(-t74, t83);
t70 = sin(t73);
t71 = cos(t73);
t64 = -t70 * t74 + t71 * t83;
t63 = 0.1e1 / t64 ^ 2;
t116 = t101 * t105;
t91 = t100 * t106 - t97 * t115;
t117 = t91 * t103;
t90 = -t100 * t104 - t97 * t114;
t77 = -t97 * t112 - t90 * t116 + t117;
t123 = t63 * t77;
t78 = t91 * t105 + (t101 * t90 + t97 * t119) * t103;
t85 = t97 * t118 - t90 * t98;
t96 = pkin(13) + qJ(5);
t94 = sin(t96);
t95 = cos(t96);
t69 = t78 * t95 + t85 * t94;
t67 = 0.1e1 / t69 ^ 2;
t68 = t78 * t94 - t85 * t95;
t122 = t67 * t68;
t82 = 0.1e1 / t83 ^ 2;
t121 = t74 * t82;
t120 = t91 * t98;
t113 = t104 * t105;
t111 = t68 ^ 2 * t67 + 0.1e1;
t110 = -t70 * t83 - t71 * t74;
t88 = (t101 * t113 + t103 * t106) * t99;
t84 = t109 * t103 + t99 * t113;
t81 = 0.1e1 / t83;
t80 = -t101 * t117 + t90 * t105;
t79 = t108 * t103 + t89 * t116;
t76 = t89 * t105 + (-t100 * t119 + t107) * t103;
t72 = 0.1e1 / (t74 ^ 2 * t82 + 0.1e1);
t66 = 0.1e1 / t69;
t65 = 0.1e1 / t111;
t62 = 0.1e1 / t64;
t61 = 0.1e1 / (t77 ^ 2 * t63 + 0.1e1);
t60 = (t88 * t121 - t79 * t81) * t72;
t59 = (t84 * t121 - t76 * t81) * t72;
t1 = [0, t60, t59, 0, 0, 0; 0 ((t90 * t103 + t91 * t116) * t62 - (t110 * t60 - t70 * t79 + t71 * t88) * t123) * t61 (t78 * t62 - (t110 * t59 - t70 * t76 + t71 * t84) * t123) * t61, 0, 0, 0; 0 ((-t95 * t120 + t80 * t94) * t66 - (t94 * t120 + t80 * t95) * t122) * t65 (t95 * t122 - t94 * t66) * t77 * t65, 0, t111 * t65, 0;];
Ja_rot  = t1;
