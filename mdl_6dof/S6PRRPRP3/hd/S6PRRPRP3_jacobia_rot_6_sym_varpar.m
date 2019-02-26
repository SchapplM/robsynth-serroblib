% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:22
% EndTime: 2019-02-26 20:02:23
% DurationCPUTime: 0.20s
% Computational Cost: add. (1132->41), mult. (1993->101), div. (87->9), fcn. (2807->13), ass. (0->58)
t100 = cos(pkin(10));
t103 = sin(qJ(2));
t101 = cos(pkin(6));
t105 = cos(qJ(2));
t108 = t101 * t105;
t98 = sin(pkin(10));
t106 = t100 * t108 - t98 * t103;
t104 = cos(qJ(3));
t102 = sin(qJ(3));
t99 = sin(pkin(6));
t113 = t102 * t99;
t109 = t101 * t103;
t90 = t100 * t109 + t98 * t105;
t85 = -t100 * t113 + t90 * t104;
t97 = pkin(11) + qJ(5);
t95 = sin(t97);
t96 = cos(t97);
t73 = t106 * t96 + t85 * t95;
t110 = t105 * t99;
t111 = t104 * t99;
t94 = t101 * t102 + t103 * t111;
t81 = t96 * t110 + t94 * t95;
t69 = atan2(-t73, t81);
t66 = sin(t69);
t67 = cos(t69);
t65 = -t66 * t73 + t67 * t81;
t64 = 0.1e1 / t65 ^ 2;
t91 = t100 * t103 + t98 * t108;
t114 = t91 * t96;
t92 = t100 * t105 - t98 * t109;
t87 = t92 * t104 + t98 * t113;
t76 = t87 * t95 - t114;
t118 = t64 * t76;
t77 = t87 * t96 + t91 * t95;
t72 = 0.1e1 / t77 ^ 2;
t86 = -t92 * t102 + t98 * t111;
t117 = t72 * t86;
t80 = 0.1e1 / t81 ^ 2;
t116 = t73 * t80;
t115 = t86 ^ 2 * t72;
t112 = t104 * t95;
t107 = -t66 * t81 - t67 * t73;
t93 = t101 * t104 - t103 * t113;
t88 = (-t103 * t96 + t105 * t112) * t99;
t84 = -t100 * t111 - t90 * t102;
t82 = -t95 * t110 + t94 * t96;
t79 = 0.1e1 / t81;
t78 = t106 * t112 - t90 * t96;
t75 = -t106 * t95 + t85 * t96;
t71 = 0.1e1 / t77;
t70 = 0.1e1 / (0.1e1 + t115);
t68 = 0.1e1 / (t73 ^ 2 * t80 + 0.1e1);
t63 = 0.1e1 / t65;
t62 = 0.1e1 / (t76 ^ 2 * t64 + 0.1e1);
t61 = (t93 * t116 - t79 * t84) * t95 * t68;
t60 = (t88 * t116 - t78 * t79) * t68;
t59 = (t82 * t116 - t75 * t79) * t68;
t1 = [0, t60, t61, 0, t59, 0; 0 ((-t91 * t112 - t92 * t96) * t63 - (t107 * t60 - t66 * t78 + t67 * t88) * t118) * t62 (t86 * t95 * t63 - ((-t66 * t84 + t67 * t93) * t95 + t107 * t61) * t118) * t62, 0 (t77 * t63 - (t107 * t59 - t66 * t75 + t67 * t82) * t118) * t62, 0; 0 (t91 * t102 * t71 - (-t104 * t114 + t92 * t95) * t117) * t70 (-t96 * t115 - t71 * t87) * t70, 0, t76 * t70 * t117, 0;];
Ja_rot  = t1;
