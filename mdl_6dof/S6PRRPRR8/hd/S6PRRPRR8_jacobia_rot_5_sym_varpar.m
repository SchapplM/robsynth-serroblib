% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR8_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:07
% EndTime: 2019-02-26 20:08:07
% DurationCPUTime: 0.18s
% Computational Cost: add. (607->40), mult. (1825->96), div. (65->9), fcn. (2498->15), ass. (0->61)
t90 = sin(pkin(7));
t91 = sin(pkin(6));
t112 = t90 * t91;
t92 = cos(pkin(12));
t104 = t92 * t112;
t93 = cos(pkin(7));
t96 = sin(qJ(3));
t110 = t93 * t96;
t100 = cos(qJ(2));
t105 = t92 * t100;
t89 = sin(pkin(12));
t94 = cos(pkin(6));
t97 = sin(qJ(2));
t83 = t94 * t105 - t89 * t97;
t106 = t89 * t100;
t109 = t94 * t97;
t84 = t92 * t109 + t106;
t99 = cos(qJ(3));
t71 = -t96 * t104 + t83 * t110 + t84 * t99;
t111 = t90 * t94;
t80 = t96 * t111 + (t100 * t110 + t97 * t99) * t91;
t69 = atan2(-t71, t80);
t66 = sin(t69);
t67 = cos(t69);
t60 = -t66 * t71 + t67 * t80;
t59 = 0.1e1 / t60 ^ 2;
t85 = -t94 * t106 - t92 * t97;
t101 = t89 * t112 + t85 * t93;
t86 = -t89 * t109 + t105;
t113 = t86 * t99;
t74 = t101 * t96 + t113;
t117 = t59 * t74;
t73 = -t101 * t99 + t86 * t96;
t81 = t89 * t91 * t93 - t85 * t90;
t95 = sin(qJ(5));
t98 = cos(qJ(5));
t65 = t73 * t95 + t81 * t98;
t63 = 0.1e1 / t65 ^ 2;
t64 = -t73 * t98 + t81 * t95;
t116 = t63 * t64;
t78 = 0.1e1 / t80 ^ 2;
t115 = t71 * t78;
t114 = t86 * t90;
t108 = t96 * t97;
t107 = t100 * t99;
t103 = t64 ^ 2 * t63 + 0.1e1;
t102 = -t66 * t80 - t67 * t71;
t82 = (-t93 * t108 + t107) * t91;
t79 = t99 * t111 + (t93 * t107 - t108) * t91;
t77 = 0.1e1 / t80;
t76 = t93 * t113 + t85 * t96;
t75 = -t84 * t110 + t83 * t99;
t70 = t84 * t96 + (-t83 * t93 + t104) * t99;
t68 = 0.1e1 / (t71 ^ 2 * t78 + 0.1e1);
t62 = 0.1e1 / t65;
t61 = 0.1e1 / t103;
t58 = 0.1e1 / t60;
t57 = 0.1e1 / (t74 ^ 2 * t59 + 0.1e1);
t56 = (t82 * t115 - t75 * t77) * t68;
t55 = (t79 * t115 + t70 * t77) * t68;
t1 = [0, t56, t55, 0, 0, 0; 0 ((-t86 * t110 + t85 * t99) * t58 - (t102 * t56 - t66 * t75 + t67 * t82) * t117) * t57 (-t73 * t58 - (t102 * t55 + t66 * t70 + t67 * t79) * t117) * t57, 0, 0, 0; 0 ((t95 * t114 - t76 * t98) * t62 - (t98 * t114 + t76 * t95) * t116) * t61 (-t95 * t116 - t98 * t62) * t74 * t61, 0, t103 * t61, 0;];
Ja_rot  = t1;
