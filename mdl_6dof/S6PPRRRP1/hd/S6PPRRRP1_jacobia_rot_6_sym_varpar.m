% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:41:33
% EndTime: 2019-02-26 19:41:34
% DurationCPUTime: 0.29s
% Computational Cost: add. (1048->43), mult. (3005->105), div. (65->9), fcn. (4119->17), ass. (0->66)
t107 = sin(pkin(11));
t112 = cos(pkin(7));
t106 = sin(pkin(12));
t110 = cos(pkin(12));
t111 = cos(pkin(11));
t113 = cos(pkin(6));
t134 = t107 * t113;
t124 = t111 * t106 + t110 * t134;
t108 = sin(pkin(7));
t109 = sin(pkin(6));
t132 = t109 * t108;
t140 = -t107 * t132 + t124 * t112;
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t129 = t111 * t113;
t125 = -t107 * t106 + t110 * t129;
t131 = t109 * t112;
t122 = -t125 * t108 - t111 * t131;
t103 = t106 * t129 + t107 * t110;
t116 = sin(qJ(3));
t119 = cos(qJ(3));
t121 = -t111 * t132 + t125 * t112;
t90 = t103 * t119 + t121 * t116;
t84 = t90 * t115 - t122 * t118;
t130 = t110 * t112;
t133 = t108 * t113;
t100 = t116 * t133 + (t106 * t119 + t116 * t130) * t109;
t102 = -t110 * t132 + t113 * t112;
t95 = t100 * t115 - t102 * t118;
t83 = atan2(-t84, t95);
t80 = sin(t83);
t81 = cos(t83);
t74 = -t80 * t84 + t81 * t95;
t73 = 0.1e1 / t74 ^ 2;
t120 = t107 * t131 + t124 * t108;
t104 = -t106 * t134 + t111 * t110;
t92 = t104 * t119 - t140 * t116;
t87 = t92 * t115 - t120 * t118;
t139 = t73 * t87;
t117 = cos(qJ(5));
t114 = sin(qJ(5));
t91 = t104 * t116 + t140 * t119;
t136 = t91 * t114;
t88 = t120 * t115 + t92 * t118;
t79 = t88 * t117 + t136;
t77 = 0.1e1 / t79 ^ 2;
t135 = t91 * t117;
t78 = t88 * t114 - t135;
t138 = t77 * t78;
t94 = 0.1e1 / t95 ^ 2;
t137 = t84 * t94;
t128 = t78 ^ 2 * t77 + 0.1e1;
t126 = -t80 * t95 - t81 * t84;
t99 = t119 * t133 + (-t106 * t116 + t119 * t130) * t109;
t96 = t100 * t118 + t102 * t115;
t93 = 0.1e1 / t95;
t89 = -t103 * t116 + t121 * t119;
t86 = t122 * t115 + t90 * t118;
t82 = 0.1e1 / (t84 ^ 2 * t94 + 0.1e1);
t76 = 0.1e1 / t79;
t75 = 0.1e1 / t128;
t72 = 0.1e1 / t74;
t71 = 0.1e1 / (t87 ^ 2 * t73 + 0.1e1);
t70 = (t99 * t137 - t89 * t93) * t82 * t115;
t69 = (t96 * t137 - t86 * t93) * t82;
t1 = [0, 0, t70, t69, 0, 0; 0, 0 (-t91 * t115 * t72 - (t126 * t70 + (-t80 * t89 + t81 * t99) * t115) * t139) * t71 (t88 * t72 - (t126 * t69 - t80 * t86 + t81 * t96) * t139) * t71, 0, 0; 0, 0 ((-t92 * t117 - t118 * t136) * t76 - (t92 * t114 - t118 * t135) * t138) * t75 (-t114 * t76 + t117 * t138) * t87 * t75, t128 * t75, 0;];
Ja_rot  = t1;
