% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:43
% EndTime: 2019-02-26 22:13:43
% DurationCPUTime: 0.32s
% Computational Cost: add. (621->32), mult. (1800->78), div. (104->9), fcn. (2548->11), ass. (0->50)
t111 = cos(qJ(3));
t107 = sin(qJ(3));
t113 = cos(qJ(1));
t121 = t113 * t107;
t109 = sin(qJ(1));
t112 = cos(qJ(2));
t123 = t109 * t112;
t100 = t111 * t123 - t121;
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t120 = t113 * t111;
t114 = t107 * t123 + t120;
t136 = -t100 * t106 + t114 * t110;
t108 = sin(qJ(2));
t115 = t106 * t111 - t107 * t110;
t98 = t115 * t108;
t82 = atan2(t136, t98);
t80 = cos(t82);
t134 = t80 * t136;
t79 = sin(t82);
t116 = -t79 * t98 + t134;
t78 = t136 * t79 + t80 * t98;
t77 = 0.1e1 / t78 ^ 2;
t101 = t109 * t107 + t112 * t120;
t117 = -t109 * t111 + t112 * t121;
t90 = t101 * t106 - t117 * t110;
t135 = t77 * t90;
t96 = 0.1e1 / t98 ^ 2;
t132 = t136 * t96;
t81 = 0.1e1 / (t136 ^ 2 * t96 + 0.1e1);
t89 = t100 * t110 + t114 * t106;
t95 = 0.1e1 / t98;
t137 = -t106 * t107 - t110 * t111;
t97 = t137 * t108;
t142 = (-t97 * t132 + t89 * t95) * t81;
t131 = t90 ^ 2 * t77;
t75 = 0.1e1 / (0.1e1 + t131);
t76 = 0.1e1 / t78;
t91 = t101 * t110 + t117 * t106;
t145 = ((t116 * t142 + t79 * t89 + t80 * t97) * t135 + t91 * t76) * t75;
t124 = t108 * t113;
t85 = 0.1e1 / t91 ^ 2;
t130 = t85 * t108 ^ 2;
t83 = 0.1e1 / (t113 ^ 2 * t130 + 0.1e1);
t138 = t90 * t83 * t85 * t124;
t99 = t115 * t112;
t92 = t109 * t98;
t84 = 0.1e1 / t91;
t74 = (-t99 * t132 + t92 * t95) * t81;
t1 = [-t90 * t95 * t81, t74, t142, 0, -t142, 0; (t136 * t76 - (-t79 + (-t95 * t134 + t79) * t81) * t131) * t75 (-(t116 * t74 + t79 * t92 + t80 * t99) * t135 - t115 * t76 * t124) * t75, -t145, 0, t145, 0; (t113 * t85 * t89 - t109 * t84) * t83 * t108 (-t113 * t130 * t137 + t112 * t84) * t83 * t113, -t138, 0, t138, 0;];
Ja_rot  = t1;
