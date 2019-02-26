% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:33
% EndTime: 2019-02-26 21:29:34
% DurationCPUTime: 0.31s
% Computational Cost: add. (1286->40), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->58)
t106 = pkin(12) + qJ(5);
t104 = sin(t106);
t105 = cos(t106);
t108 = sin(pkin(6));
t116 = cos(qJ(1));
t122 = t108 * t116;
t113 = sin(qJ(1));
t110 = cos(pkin(6));
t107 = sin(pkin(11));
t109 = cos(pkin(11));
t112 = sin(qJ(2));
t115 = cos(qJ(2));
t118 = t115 * t107 + t112 * t109;
t98 = t118 * t110;
t99 = t112 * t107 - t115 * t109;
t87 = -t113 * t99 + t116 * t98;
t80 = t87 * t104 + t105 * t122;
t97 = t118 * t108;
t93 = t97 * t104 - t110 * t105;
t79 = atan2(-t80, t93);
t76 = sin(t79);
t77 = cos(t79);
t70 = -t76 * t80 + t77 * t93;
t69 = 0.1e1 / t70 ^ 2;
t119 = -t113 * t98 - t116 * t99;
t123 = t108 * t113;
t84 = t104 * t119 - t105 * t123;
t130 = t69 * t84;
t114 = cos(qJ(6));
t111 = sin(qJ(6));
t117 = t99 * t110;
t89 = t113 * t117 - t116 * t118;
t125 = t89 * t111;
t85 = t104 * t123 + t105 * t119;
t75 = t85 * t114 - t125;
t73 = 0.1e1 / t75 ^ 2;
t124 = t89 * t114;
t74 = t85 * t111 + t124;
t129 = t73 * t74;
t128 = t77 * t80;
t92 = 0.1e1 / t93 ^ 2;
t127 = t80 * t92;
t126 = t84 ^ 2 * t69;
t121 = t74 ^ 2 * t73 + 0.1e1;
t82 = -t104 * t122 + t87 * t105;
t120 = -t76 * t93 - t128;
t96 = t99 * t108;
t94 = t110 * t104 + t97 * t105;
t91 = 0.1e1 / t93;
t86 = -t113 * t118 - t116 * t117;
t78 = 0.1e1 / (t80 ^ 2 * t92 + 0.1e1);
t72 = 0.1e1 / t75;
t71 = 0.1e1 / t121;
t68 = 0.1e1 / t70;
t67 = 0.1e1 / (0.1e1 + t126);
t66 = (-t96 * t127 - t86 * t91) * t78 * t104;
t65 = (t94 * t127 - t82 * t91) * t78;
t1 = [-t84 * t91 * t78, t66, 0, 0, t65, 0; (-t80 * t68 - (-t76 + (t91 * t128 + t76) * t78) * t126) * t67 (t89 * t104 * t68 - (t120 * t66 + (-t76 * t86 - t77 * t96) * t104) * t130) * t67, 0, 0 (t85 * t68 - (t120 * t65 - t76 * t82 + t77 * t94) * t130) * t67, 0; ((-t111 * t82 - t86 * t114) * t72 - (t86 * t111 - t114 * t82) * t129) * t71 ((t105 * t125 - t114 * t119) * t72 - (t105 * t124 + t111 * t119) * t129) * t71, 0, 0 (-t111 * t72 + t114 * t129) * t84 * t71, t121 * t71;];
Ja_rot  = t1;
