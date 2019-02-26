% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:56
% EndTime: 2019-02-26 20:12:56
% DurationCPUTime: 0.22s
% Computational Cost: add. (655->41), mult. (1825->97), div. (65->9), fcn. (2498->15), ass. (0->61)
t101 = cos(pkin(12));
t102 = cos(pkin(7));
t105 = sin(qJ(2));
t103 = cos(pkin(6));
t107 = cos(qJ(2));
t116 = t103 * t107;
t98 = sin(pkin(12));
t109 = t101 * t116 - t105 * t98;
t100 = sin(pkin(6));
t99 = sin(pkin(7));
t121 = t100 * t99;
t126 = -t101 * t121 + t102 * t109;
t104 = sin(qJ(3));
t106 = cos(qJ(3));
t117 = t103 * t105;
t90 = t101 * t117 + t107 * t98;
t75 = t90 * t104 - t106 * t126;
t118 = t102 * t106;
t120 = t103 * t99;
t84 = -t106 * t120 + (t104 * t105 - t107 * t118) * t100;
t74 = atan2(-t75, t84);
t71 = sin(t74);
t72 = cos(t74);
t65 = -t71 * t75 + t72 * t84;
t64 = 0.1e1 / t65 ^ 2;
t113 = t98 * t121;
t92 = t101 * t107 - t117 * t98;
t119 = t92 * t104;
t91 = -t101 * t105 - t116 * t98;
t78 = -t106 * t113 - t118 * t91 + t119;
t125 = t64 * t78;
t79 = t92 * t106 + (t102 * t91 + t113) * t104;
t86 = t100 * t102 * t98 - t91 * t99;
t97 = qJ(4) + pkin(13);
t95 = sin(t97);
t96 = cos(t97);
t70 = t79 * t96 + t86 * t95;
t68 = 0.1e1 / t70 ^ 2;
t69 = t79 * t95 - t86 * t96;
t124 = t68 * t69;
t83 = 0.1e1 / t84 ^ 2;
t123 = t75 * t83;
t122 = t92 * t99;
t115 = t104 * t107;
t114 = t105 * t106;
t111 = t68 * t69 ^ 2 + 0.1e1;
t110 = -t71 * t84 - t72 * t75;
t89 = (t102 * t114 + t115) * t100;
t85 = t104 * t120 + (t102 * t115 + t114) * t100;
t82 = 0.1e1 / t84;
t81 = -t102 * t119 + t106 * t91;
t80 = t104 * t109 + t118 * t90;
t77 = t104 * t126 + t90 * t106;
t73 = 0.1e1 / (t75 ^ 2 * t83 + 0.1e1);
t67 = 0.1e1 / t70;
t66 = 0.1e1 / t111;
t63 = 0.1e1 / t65;
t62 = 0.1e1 / (t64 * t78 ^ 2 + 0.1e1);
t61 = (t123 * t89 - t80 * t82) * t73;
t60 = (t123 * t85 - t77 * t82) * t73;
t1 = [0, t61, t60, 0, 0, 0; 0 ((t104 * t91 + t118 * t92) * t63 - (t110 * t61 - t71 * t80 + t72 * t89) * t125) * t62 (t79 * t63 - (t110 * t60 - t71 * t77 + t72 * t85) * t125) * t62, 0, 0, 0; 0 ((-t122 * t96 + t81 * t95) * t67 - (t122 * t95 + t81 * t96) * t124) * t66 (t124 * t96 - t67 * t95) * t78 * t66, t111 * t66, 0, 0;];
Ja_rot  = t1;
