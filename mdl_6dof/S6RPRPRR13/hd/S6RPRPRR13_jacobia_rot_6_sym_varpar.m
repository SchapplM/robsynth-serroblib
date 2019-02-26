% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR13_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:43
% EndTime: 2019-02-26 20:55:43
% DurationCPUTime: 0.45s
% Computational Cost: add. (1440->49), mult. (4121->118), div. (85->9), fcn. (5660->17), ass. (0->69)
t112 = cos(pkin(12));
t109 = sin(pkin(12));
t118 = sin(qJ(1));
t133 = t118 * t109;
t114 = cos(pkin(6));
t122 = cos(qJ(1));
t134 = t114 * t122;
t106 = -t112 * t134 + t133;
t110 = sin(pkin(7));
t113 = cos(pkin(7));
t111 = sin(pkin(6));
t136 = t111 * t122;
t101 = -t106 * t110 + t113 * t136;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t132 = t118 * t112;
t107 = t109 * t134 + t132;
t117 = sin(qJ(3));
t121 = cos(qJ(3));
t128 = t106 * t113 + t110 * t136;
t123 = t107 * t117 + t128 * t121;
t86 = -t101 * t120 + t123 * t116;
t151 = t101 * t116 + t123 * t120;
t127 = -t109 * t122 - t114 * t132;
t137 = t111 * t118;
t148 = t110 * t137 + t127 * t113;
t97 = -t107 * t121 + t128 * t117;
t105 = -t110 * t111 * t112 + t113 * t114;
t135 = t112 * t113;
t138 = t110 * t114;
t99 = -t121 * t138 + (t109 * t117 - t121 * t135) * t111;
t93 = t105 * t116 - t120 * t99;
t82 = atan2(t151, t93);
t79 = sin(t82);
t80 = cos(t82);
t73 = t151 * t79 + t80 * t93;
t72 = 0.1e1 / t73 ^ 2;
t103 = -t127 * t110 + t113 * t137;
t126 = t112 * t122 - t114 * t133;
t124 = t126 * t117 - t148 * t121;
t87 = t103 * t116 - t124 * t120;
t147 = t72 * t87;
t146 = t72 * t87 ^ 2;
t119 = cos(qJ(6));
t115 = sin(qJ(6));
t98 = t148 * t117 + t126 * t121;
t142 = t115 * t98;
t88 = t103 * t120 + t124 * t116;
t78 = t119 * t88 + t142;
t76 = 0.1e1 / t78 ^ 2;
t141 = t119 * t98;
t77 = t115 * t88 - t141;
t145 = t76 * t77;
t144 = t80 * t151;
t92 = 0.1e1 / t93 ^ 2;
t143 = t151 * t92;
t131 = t76 * t77 ^ 2 + 0.1e1;
t129 = -t79 * t93 + t144;
t100 = t117 * t138 + (t109 * t121 + t117 * t135) * t111;
t94 = t105 * t120 + t116 * t99;
t91 = 0.1e1 / t93;
t81 = 0.1e1 / (t151 ^ 2 * t92 + 0.1e1);
t75 = 0.1e1 / t78;
t74 = 0.1e1 / t131;
t71 = 0.1e1 / t73;
t70 = 0.1e1 / (0.1e1 + t146);
t69 = (t100 * t143 - t91 * t97) * t81 * t120;
t68 = (-t94 * t143 - t86 * t91) * t81;
t1 = [-t87 * t91 * t81, 0, t69, 0, t68, 0; (t151 * t71 - (-t79 + (-t91 * t144 + t79) * t81) * t146) * t70, 0 (-t98 * t120 * t71 - (t129 * t69 + (-t100 * t80 - t79 * t97) * t120) * t147) * t70, 0 (t88 * t71 - (t129 * t68 - t79 * t86 + t80 * t94) * t147) * t70, 0; ((-t115 * t86 - t119 * t97) * t75 - (t115 * t97 - t119 * t86) * t145) * t74, 0 ((t116 * t142 + t124 * t119) * t75 - (-t124 * t115 + t116 * t141) * t145) * t74, 0 (-t115 * t75 + t119 * t145) * t87 * t74, t131 * t74;];
Ja_rot  = t1;
