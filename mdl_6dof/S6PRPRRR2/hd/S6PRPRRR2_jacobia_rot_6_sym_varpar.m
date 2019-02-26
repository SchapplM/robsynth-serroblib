% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR2
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
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:12
% EndTime: 2019-02-26 19:54:12
% DurationCPUTime: 0.17s
% Computational Cost: add. (753->33), mult. (1880->84), div. (70->9), fcn. (2635->15), ass. (0->56)
t105 = cos(pkin(11));
t107 = sin(qJ(4));
t103 = sin(pkin(6));
t109 = cos(qJ(4));
t116 = t103 * t109;
t102 = sin(pkin(11));
t106 = cos(pkin(6));
t101 = sin(pkin(12));
t104 = cos(pkin(12));
t108 = sin(qJ(2));
t110 = cos(qJ(2));
t112 = t110 * t101 + t108 * t104;
t93 = t112 * t106;
t94 = t108 * t101 - t110 * t104;
t82 = -t102 * t94 + t105 * t93;
t76 = t105 * t116 + t82 * t107;
t92 = t112 * t103;
t88 = -t106 * t109 + t92 * t107;
t75 = atan2(-t76, t88);
t72 = sin(t75);
t73 = cos(t75);
t66 = -t72 * t76 + t73 * t88;
t65 = 0.1e1 / t66 ^ 2;
t113 = -t102 * t93 - t105 * t94;
t79 = -t102 * t116 + t107 * t113;
t121 = t65 * t79;
t117 = t103 * t107;
t80 = t102 * t117 + t109 * t113;
t111 = t94 * t106;
t84 = t102 * t111 - t105 * t112;
t100 = qJ(5) + qJ(6);
t98 = sin(t100);
t99 = cos(t100);
t71 = t80 * t99 - t84 * t98;
t69 = 0.1e1 / t71 ^ 2;
t70 = t80 * t98 + t84 * t99;
t120 = t69 * t70;
t87 = 0.1e1 / t88 ^ 2;
t119 = t76 * t87;
t118 = t109 * t84;
t115 = t70 ^ 2 * t69 + 0.1e1;
t114 = -t72 * t88 - t73 * t76;
t91 = t94 * t103;
t89 = t106 * t107 + t92 * t109;
t86 = 0.1e1 / t88;
t81 = -t102 * t112 - t105 * t111;
t78 = -t105 * t117 + t82 * t109;
t74 = 0.1e1 / (t76 ^ 2 * t87 + 0.1e1);
t68 = 0.1e1 / t71;
t67 = 0.1e1 / t115;
t64 = 0.1e1 / t66;
t63 = 0.1e1 / (t79 ^ 2 * t65 + 0.1e1);
t62 = (-t91 * t119 - t81 * t86) * t74 * t107;
t61 = (t89 * t119 - t78 * t86) * t74;
t60 = t115 * t67;
t1 = [0, t62, 0, t61, 0, 0; 0 (t84 * t107 * t64 - (t114 * t62 + (-t72 * t81 - t73 * t91) * t107) * t121) * t63, 0 (t80 * t64 - (t114 * t61 - t72 * t78 + t73 * t89) * t121) * t63, 0, 0; 0 ((-t113 * t99 + t98 * t118) * t68 - (t113 * t98 + t99 * t118) * t120) * t67, 0 (t99 * t120 - t68 * t98) * t79 * t67, t60, t60;];
Ja_rot  = t1;
