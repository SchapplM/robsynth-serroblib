% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR15_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:20
% EndTime: 2019-02-26 22:24:20
% DurationCPUTime: 0.24s
% Computational Cost: add. (662->45), mult. (2020->99), div. (77->9), fcn. (2784->13), ass. (0->59)
t90 = sin(pkin(6));
t98 = cos(qJ(1));
t111 = t90 * t98;
t89 = sin(pkin(7));
t101 = t89 * t111;
t91 = cos(pkin(7));
t96 = cos(qJ(3));
t109 = t91 * t96;
t94 = sin(qJ(2));
t104 = t98 * t94;
t95 = sin(qJ(1));
t97 = cos(qJ(2));
t105 = t95 * t97;
t92 = cos(pkin(6));
t83 = t104 * t92 + t105;
t93 = sin(qJ(3));
t114 = t83 * t93;
t103 = t98 * t97;
t106 = t95 * t94;
t82 = -t103 * t92 + t106;
t66 = t101 * t96 + t109 * t82 + t114;
t113 = t89 * t92;
t74 = -t96 * t113 + (-t109 * t97 + t93 * t94) * t90;
t64 = atan2(-t66, t74);
t61 = sin(t64);
t62 = cos(t64);
t60 = -t61 * t66 + t62 * t74;
t59 = 0.1e1 / t60 ^ 2;
t112 = t90 * t95;
t102 = t89 * t112;
t84 = -t105 * t92 - t104;
t85 = -t106 * t92 + t103;
t69 = -t102 * t96 - t109 * t84 + t85 * t93;
t119 = t59 * t69;
t118 = t62 * t66;
t73 = 0.1e1 / t74 ^ 2;
t117 = t66 * t73;
t116 = t69 ^ 2 * t59;
t70 = t85 * t96 + (t84 * t91 + t102) * t93;
t78 = t112 * t91 - t84 * t89;
t77 = 0.1e1 / t78 ^ 2;
t115 = t70 * t77;
t110 = t91 * t93;
t108 = t93 * t97;
t107 = t94 * t96;
t100 = -t61 * t74 - t118;
t99 = t101 * t93 + t110 * t82 - t83 * t96;
t81 = (t107 * t91 + t108) * t90;
t76 = 0.1e1 / t78;
t75 = t93 * t113 + (t108 * t91 + t107) * t90;
t72 = 0.1e1 / t74;
t71 = t109 * t83 - t82 * t93;
t65 = 0.1e1 / (t70 ^ 2 * t77 + 0.1e1);
t63 = 0.1e1 / (t66 ^ 2 * t73 + 0.1e1);
t58 = 0.1e1 / t60;
t57 = 0.1e1 / (0.1e1 + t116);
t56 = (t117 * t81 - t71 * t72) * t63;
t55 = (t117 * t75 + t72 * t99) * t63;
t1 = [-t69 * t72 * t63, t56, t55, 0, 0, 0; ((-t114 + (-t82 * t91 - t101) * t96) * t58 - (-t61 + (t118 * t72 + t61) * t63) * t116) * t57 ((t109 * t85 + t84 * t93) * t58 - (t100 * t56 - t61 * t71 + t62 * t81) * t119) * t57 (t70 * t58 - (t100 * t55 + t61 * t99 + t62 * t75) * t119) * t57, 0, 0, 0; (t99 * t76 - (t111 * t91 - t82 * t89) * t115) * t65 ((-t110 * t85 + t84 * t96) * t76 - t85 * t89 * t115) * t65, -t69 * t76 * t65, 0, 0, 0;];
Ja_rot  = t1;
