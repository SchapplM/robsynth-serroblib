% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14V3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobia_rot_6_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:09
% EndTime: 2019-04-12 15:12:09
% DurationCPUTime: 0.24s
% Computational Cost: add. (609->45), mult. (1737->115), div. (115->9), fcn. (2479->13), ass. (0->58)
t95 = sin(qJ(2));
t98 = cos(qJ(5));
t110 = t95 * t98;
t101 = cos(qJ(1));
t100 = cos(qJ(2));
t99 = cos(qJ(4));
t107 = t100 * t99;
t94 = sin(qJ(4));
t96 = sin(qJ(1));
t84 = -t101 * t94 + t96 * t107;
t93 = sin(qJ(5));
t71 = -t96 * t110 + t84 * t93;
t109 = t95 * t99;
t81 = t100 * t98 + t93 * t109;
t70 = atan2(-t71, t81);
t67 = sin(t70);
t68 = cos(t70);
t61 = -t67 * t71 + t68 * t81;
t60 = 0.1e1 / t61 ^ 2;
t106 = t101 * t95;
t105 = t100 * t101;
t108 = t96 * t94;
t87 = t99 * t105 + t108;
t75 = -t98 * t106 + t87 * t93;
t117 = t60 * t75;
t76 = t93 * t106 + t87 * t98;
t86 = t94 * t105 - t96 * t99;
t92 = sin(qJ(6));
t97 = cos(qJ(6));
t66 = t76 * t97 + t86 * t92;
t64 = 0.1e1 / t66 ^ 2;
t65 = t76 * t92 - t86 * t97;
t116 = t64 * t65;
t115 = t68 * t71;
t80 = 0.1e1 / t81 ^ 2;
t114 = t71 * t80;
t113 = t75 ^ 2 * t60;
t112 = t86 * t98;
t111 = t94 * t95;
t104 = t94 * t106;
t103 = t65 ^ 2 * t64 + 0.1e1;
t102 = -t67 * t81 - t115;
t73 = t96 * t95 * t93 + t84 * t98;
t82 = -t100 * t93 + t98 * t109;
t85 = t93 * t107 - t110;
t83 = -t100 * t108 - t101 * t99;
t79 = 0.1e1 / t81;
t78 = t82 * t101;
t77 = t81 * t96;
t69 = 0.1e1 / (t71 ^ 2 * t80 + 0.1e1);
t63 = 0.1e1 / t66;
t62 = 0.1e1 / t103;
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (0.1e1 + t113);
t57 = (-t111 * t114 - t79 * t83) * t93 * t69;
t56 = (t85 * t114 + t77 * t79) * t69;
t55 = (t82 * t114 - t73 * t79) * t69;
t1 = [-t75 * t79 * t69, t56, 0, t57, t55, 0; (-t71 * t59 - (-t67 + (t79 * t115 + t67) * t69) * t113) * t58 (-(t102 * t56 + t67 * t77 + t68 * t85) * t117 - t81 * t59 * t101) * t58, 0 (-t86 * t93 * t59 - ((-t68 * t111 - t67 * t83) * t93 + t102 * t57) * t117) * t58 (t76 * t59 - (t102 * t55 - t67 * t73 + t68 * t82) * t117) * t58, 0; ((-t73 * t92 - t83 * t97) * t63 - (-t73 * t97 + t83 * t92) * t116) * t62 ((t97 * t104 - t78 * t92) * t63 - (-t92 * t104 - t78 * t97) * t116) * t62, 0 ((-t92 * t112 - t87 * t97) * t63 - (-t97 * t112 + t87 * t92) * t116) * t62 (t97 * t116 - t92 * t63) * t75 * t62, t103 * t62;];
Ja_rot  = t1;
