% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:31
% EndTime: 2019-02-26 20:03:31
% DurationCPUTime: 0.18s
% Computational Cost: add. (716->41), mult. (1999->102), div. (87->9), fcn. (2816->13), ass. (0->57)
t93 = cos(pkin(6));
t96 = sin(qJ(2));
t106 = t93 * t96;
t90 = sin(pkin(10));
t92 = cos(pkin(10));
t99 = cos(qJ(2));
t101 = t92 * t106 + t90 * t99;
t91 = sin(pkin(6));
t98 = cos(qJ(3));
t107 = t91 * t98;
t95 = sin(qJ(3));
t100 = t101 * t95 + t92 * t107;
t105 = t93 * t99;
t85 = -t92 * t105 + t90 * t96;
t94 = sin(qJ(5));
t97 = cos(qJ(5));
t69 = -t100 * t97 + t85 * t94;
t108 = t91 * t95;
t88 = t96 * t108 - t93 * t98;
t82 = -t91 * t99 * t94 - t88 * t97;
t66 = atan2(-t69, t82);
t63 = sin(t66);
t64 = cos(t66);
t61 = -t63 * t69 + t64 * t82;
t60 = 0.1e1 / t61 ^ 2;
t86 = t90 * t105 + t92 * t96;
t109 = t86 * t94;
t87 = -t90 * t106 + t92 * t99;
t77 = -t90 * t107 + t87 * t95;
t72 = -t77 * t97 + t109;
t113 = t60 * t72;
t73 = t77 * t94 + t86 * t97;
t68 = 0.1e1 / t73 ^ 2;
t78 = t90 * t108 + t87 * t98;
t112 = t68 * t78;
t81 = 0.1e1 / t82 ^ 2;
t111 = t69 * t81;
t110 = t78 ^ 2 * t68;
t104 = t95 * t97;
t103 = t97 * t99;
t102 = -t63 * t82 - t64 * t69;
t89 = t96 * t107 + t93 * t95;
t84 = (-t95 * t103 + t94 * t96) * t91;
t83 = -t91 * t103 + t88 * t94;
t80 = 0.1e1 / t82;
t76 = t101 * t98 - t92 * t108;
t74 = t101 * t94 + t85 * t104;
t71 = t100 * t94 + t85 * t97;
t67 = 0.1e1 / t73;
t65 = 0.1e1 / (t69 ^ 2 * t81 + 0.1e1);
t62 = 0.1e1 / (0.1e1 + t110);
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (t72 ^ 2 * t60 + 0.1e1);
t57 = (-t89 * t111 + t76 * t80) * t97 * t65;
t56 = (t84 * t111 - t74 * t80) * t65;
t55 = (t83 * t111 - t71 * t80) * t65;
t1 = [0, t56, t57, 0, t55, 0; 0 ((t86 * t104 + t87 * t94) * t59 - (t102 * t56 - t63 * t74 + t64 * t84) * t113) * t58 (-t78 * t97 * t59 - ((t63 * t76 - t64 * t89) * t97 + t102 * t57) * t113) * t58, 0 (t73 * t59 - (t102 * t55 - t63 * t71 + t64 * t83) * t113) * t58, 0; 0 (t86 * t98 * t67 + (-t95 * t109 + t87 * t97) * t112) * t62 (t94 * t110 + t67 * t77) * t62, 0, -t72 * t62 * t112, 0;];
Ja_rot  = t1;
