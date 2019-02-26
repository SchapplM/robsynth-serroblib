% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:52
% DurationCPUTime: 0.39s
% Computational Cost: add. (1084->41), mult. (2952->94), div. (95->9), fcn. (4136->15), ass. (0->60)
t113 = cos(pkin(11));
t121 = cos(qJ(3));
t111 = sin(pkin(11));
t122 = cos(qJ(2));
t114 = cos(pkin(6));
t118 = sin(qJ(2));
t130 = t114 * t118;
t124 = t111 * t122 + t113 * t130;
t112 = sin(pkin(6));
t117 = sin(qJ(3));
t132 = t112 * t117;
t100 = -t113 * t132 + t124 * t121;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t131 = t112 * t121;
t123 = t113 * t131 + t124 * t117;
t85 = t100 * t116 - t123 * t120;
t107 = -t114 * t121 + t118 * t132;
t108 = t114 * t117 + t118 * t131;
t96 = -t107 * t120 + t108 * t116;
t79 = atan2(-t85, t96);
t76 = sin(t79);
t77 = cos(t79);
t127 = -t76 * t96 - t77 * t85;
t74 = -t76 * t85 + t77 * t96;
t73 = 0.1e1 / t74 ^ 2;
t106 = -t111 * t130 + t113 * t122;
t101 = t106 * t117 - t111 * t131;
t102 = t106 * t121 + t111 * t132;
t89 = -t101 * t120 + t102 * t116;
t137 = t73 * t89;
t94 = 0.1e1 / t96 ^ 2;
t135 = t85 * t94;
t78 = 0.1e1 / (t85 ^ 2 * t94 + 0.1e1);
t87 = t100 * t120 + t123 * t116;
t93 = 0.1e1 / t96;
t97 = t107 * t116 + t108 * t120;
t143 = (-t97 * t135 + t87 * t93) * t78;
t71 = 0.1e1 / (t89 ^ 2 * t73 + 0.1e1);
t72 = 0.1e1 / t74;
t90 = t101 * t116 + t102 * t120;
t146 = ((t127 * t143 + t76 * t87 - t77 * t97) * t137 + t90 * t72) * t71;
t115 = sin(qJ(6));
t119 = cos(qJ(6));
t129 = t114 * t122;
t105 = -t111 * t129 - t113 * t118;
t83 = t105 * t115 + t90 * t119;
t81 = 0.1e1 / t83 ^ 2;
t82 = -t105 * t119 + t90 * t115;
t136 = t81 * t82;
t128 = t82 ^ 2 * t81 + 0.1e1;
t75 = 0.1e1 / t128;
t80 = 0.1e1 / t83;
t138 = (-t115 * t80 + t119 * t136) * t89 * t75;
t126 = t116 * t121 - t117 * t120;
t103 = t126 * t122 * t112;
t92 = (t116 * t117 + t120 * t121) * t105;
t91 = t126 * (-t111 * t118 + t113 * t129);
t70 = (t103 * t135 - t91 * t93) * t78;
t1 = [0, t70, t143, 0, -t143, 0; 0 (-(t77 * t103 + t127 * t70 - t76 * t91) * t137 + t126 * t72 * t105) * t71, -t146, 0, t146, 0; 0 ((t106 * t119 + t92 * t115) * t80 - (-t106 * t115 + t92 * t119) * t136) * t75, -t138, 0, t138, t128 * t75;];
Ja_rot  = t1;
