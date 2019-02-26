% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:43
% EndTime: 2019-02-26 19:53:43
% DurationCPUTime: 0.24s
% Computational Cost: add. (1379->33), mult. (2514->83), div. (95->9), fcn. (3534->15), ass. (0->58)
t109 = qJ(4) + qJ(5);
t107 = sin(t109);
t108 = cos(t109);
t112 = sin(pkin(6));
t114 = cos(pkin(11));
t125 = t112 * t114;
t115 = cos(pkin(6));
t110 = sin(pkin(12));
t113 = cos(pkin(12));
t117 = sin(qJ(2));
t119 = cos(qJ(2));
t121 = t119 * t110 + t117 * t113;
t102 = t121 * t115;
t103 = t117 * t110 - t119 * t113;
t111 = sin(pkin(11));
t91 = t114 * t102 - t111 * t103;
t85 = t91 * t107 + t108 * t125;
t101 = t121 * t112;
t97 = t101 * t107 - t115 * t108;
t84 = atan2(-t85, t97);
t81 = sin(t84);
t82 = cos(t84);
t75 = -t81 * t85 + t82 * t97;
t74 = 0.1e1 / t75 ^ 2;
t122 = -t111 * t102 - t114 * t103;
t126 = t111 * t112;
t88 = t107 * t122 - t108 * t126;
t131 = t74 * t88;
t118 = cos(qJ(6));
t116 = sin(qJ(6));
t120 = t103 * t115;
t93 = t111 * t120 - t114 * t121;
t128 = t93 * t116;
t89 = t107 * t126 + t108 * t122;
t80 = t89 * t118 - t128;
t78 = 0.1e1 / t80 ^ 2;
t127 = t93 * t118;
t79 = t89 * t116 + t127;
t130 = t78 * t79;
t96 = 0.1e1 / t97 ^ 2;
t129 = t85 * t96;
t124 = t79 ^ 2 * t78 + 0.1e1;
t123 = -t81 * t97 - t82 * t85;
t100 = t103 * t112;
t98 = t101 * t108 + t115 * t107;
t95 = 0.1e1 / t97;
t90 = -t111 * t121 - t114 * t120;
t87 = -t107 * t125 + t91 * t108;
t83 = 0.1e1 / (t85 ^ 2 * t96 + 0.1e1);
t77 = 0.1e1 / t80;
t76 = 0.1e1 / t124;
t73 = 0.1e1 / t75;
t72 = 0.1e1 / (t88 ^ 2 * t74 + 0.1e1);
t71 = (-t100 * t129 - t90 * t95) * t83 * t107;
t70 = (t98 * t129 - t87 * t95) * t83;
t69 = (-t116 * t77 + t118 * t130) * t88 * t76;
t68 = (t89 * t73 - (t123 * t70 - t81 * t87 + t82 * t98) * t131) * t72;
t1 = [0, t71, 0, t70, t70, 0; 0 (t93 * t107 * t73 - (t123 * t71 + (-t100 * t82 - t81 * t90) * t107) * t131) * t72, 0, t68, t68, 0; 0 ((t108 * t128 - t118 * t122) * t77 - (t108 * t127 + t116 * t122) * t130) * t76, 0, t69, t69, t124 * t76;];
Ja_rot  = t1;
