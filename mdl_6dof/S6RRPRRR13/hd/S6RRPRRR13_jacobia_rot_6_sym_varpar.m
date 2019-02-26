% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR13_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:01:09
% EndTime: 2019-02-26 22:01:09
% DurationCPUTime: 0.18s
% Computational Cost: add. (570->37), mult. (1383->91), div. (90->9), fcn. (1970->13), ass. (0->56)
t89 = sin(pkin(6));
t96 = cos(qJ(1));
t105 = t89 * t96;
t95 = cos(qJ(2));
t101 = t96 * t95;
t92 = sin(qJ(2));
t93 = sin(qJ(1));
t104 = t93 * t92;
t90 = cos(pkin(6));
t82 = -t90 * t101 + t104;
t91 = sin(qJ(4));
t94 = cos(qJ(4));
t74 = t91 * t105 + t82 * t94;
t106 = t89 * t95;
t80 = t94 * t106 + t90 * t91;
t71 = atan2(t74, t80);
t68 = sin(t71);
t69 = cos(t71);
t62 = t68 * t74 + t69 * t80;
t61 = 0.1e1 / t62 ^ 2;
t107 = t89 * t93;
t102 = t96 * t92;
t103 = t93 * t95;
t97 = t90 * t103 + t102;
t72 = t91 * t107 - t94 * t97;
t114 = t61 * t72;
t73 = t94 * t107 + t91 * t97;
t84 = -t90 * t104 + t101;
t88 = qJ(5) + qJ(6);
t86 = sin(t88);
t87 = cos(t88);
t67 = t73 * t87 + t84 * t86;
t65 = 0.1e1 / t67 ^ 2;
t66 = t73 * t86 - t84 * t87;
t113 = t65 * t66;
t112 = t69 * t74;
t111 = t72 ^ 2 * t61;
t79 = 0.1e1 / t80 ^ 2;
t110 = t74 * t79;
t109 = t84 * t91;
t108 = t89 * t92;
t100 = t66 ^ 2 * t65 + 0.1e1;
t99 = -t68 * t80 + t112;
t98 = t94 * t105 - t82 * t91;
t83 = t90 * t102 + t103;
t81 = -t91 * t106 + t90 * t94;
t78 = 0.1e1 / t80;
t70 = 0.1e1 / (t74 ^ 2 * t79 + 0.1e1);
t64 = 0.1e1 / t67;
t63 = 0.1e1 / t100;
t60 = 0.1e1 / t62;
t59 = 0.1e1 / (0.1e1 + t111);
t58 = (t108 * t110 + t78 * t83) * t94 * t70;
t57 = (-t81 * t110 + t78 * t98) * t70;
t56 = t100 * t63;
t1 = [-t72 * t78 * t70, t58, 0, t57, 0, 0; (t74 * t60 - (-t68 + (-t78 * t112 + t68) * t70) * t111) * t59 (-t84 * t94 * t60 - ((-t69 * t108 + t68 * t83) * t94 + t99 * t58) * t114) * t59, 0 (t73 * t60 - (t99 * t57 + t68 * t98 + t69 * t81) * t114) * t59, 0, 0; ((t83 * t87 + t86 * t98) * t64 - (-t83 * t86 + t87 * t98) * t113) * t63 ((t86 * t109 + t87 * t97) * t64 - (t87 * t109 - t86 * t97) * t113) * t63, 0 (t87 * t113 - t86 * t64) * t72 * t63, t56, t56;];
Ja_rot  = t1;
