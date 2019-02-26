% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:14
% EndTime: 2019-02-26 20:19:15
% DurationCPUTime: 0.17s
% Computational Cost: add. (1071->32), mult. (1470->80), div. (100->9), fcn. (2088->13), ass. (0->55)
t110 = sin(pkin(12));
t112 = cos(pkin(12));
t115 = cos(qJ(2));
t113 = cos(pkin(6));
t114 = sin(qJ(2));
t119 = t113 * t114;
t100 = t110 * t115 + t112 * t119;
t109 = qJ(3) + qJ(4);
t105 = sin(t109);
t107 = cos(t109);
t111 = sin(pkin(6));
t122 = t111 * t112;
t90 = t100 * t105 + t107 * t122;
t121 = t111 * t114;
t97 = t105 * t121 - t113 * t107;
t89 = atan2(-t90, t97);
t86 = sin(t89);
t87 = cos(t89);
t80 = -t86 * t90 + t87 * t97;
t79 = 0.1e1 / t80 ^ 2;
t102 = -t110 * t119 + t112 * t115;
t123 = t110 * t111;
t93 = t102 * t105 - t107 * t123;
t127 = t79 * t93;
t118 = t113 * t115;
t101 = t110 * t118 + t112 * t114;
t108 = qJ(5) + qJ(6);
t104 = sin(t108);
t106 = cos(t108);
t94 = t102 * t107 + t105 * t123;
t85 = t101 * t104 + t94 * t106;
t83 = 0.1e1 / t85 ^ 2;
t84 = -t101 * t106 + t94 * t104;
t126 = t83 * t84;
t96 = 0.1e1 / t97 ^ 2;
t125 = t90 * t96;
t124 = t101 * t107;
t120 = t111 * t115;
t117 = t84 ^ 2 * t83 + 0.1e1;
t116 = -t86 * t97 - t87 * t90;
t99 = -t110 * t114 + t112 * t118;
t98 = t113 * t105 + t107 * t121;
t95 = 0.1e1 / t97;
t92 = t100 * t107 - t105 * t122;
t88 = 0.1e1 / (t90 ^ 2 * t96 + 0.1e1);
t82 = 0.1e1 / t85;
t81 = 0.1e1 / t117;
t78 = 0.1e1 / t80;
t77 = 0.1e1 / (t93 ^ 2 * t79 + 0.1e1);
t76 = (t120 * t125 - t95 * t99) * t88 * t105;
t75 = (t98 * t125 - t92 * t95) * t88;
t74 = t117 * t81;
t73 = (-t104 * t82 + t106 * t126) * t93 * t81;
t72 = (t94 * t78 - (t116 * t75 - t86 * t92 + t87 * t98) * t127) * t77;
t1 = [0, t76, t75, t75, 0, 0; 0 (-t101 * t105 * t78 - (t116 * t76 + (t87 * t120 - t86 * t99) * t105) * t127) * t77, t72, t72, 0, 0; 0 ((-t102 * t106 - t104 * t124) * t82 - (t102 * t104 - t106 * t124) * t126) * t81, t73, t73, t74, t74;];
Ja_rot  = t1;
