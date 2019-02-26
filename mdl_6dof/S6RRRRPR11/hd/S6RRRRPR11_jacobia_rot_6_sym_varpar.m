% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:21
% EndTime: 2019-02-26 22:36:21
% DurationCPUTime: 0.23s
% Computational Cost: add. (650->38), mult. (1383->91), div. (90->9), fcn. (1970->13), ass. (0->56)
t100 = cos(qJ(1));
t93 = sin(pkin(6));
t105 = t100 * t93;
t96 = sin(qJ(2));
t104 = t100 * t96;
t97 = sin(qJ(1));
t99 = cos(qJ(2));
t106 = t97 * t99;
t94 = cos(pkin(6));
t85 = t94 * t104 + t106;
t95 = sin(qJ(3));
t98 = cos(qJ(3));
t74 = t98 * t105 + t85 * t95;
t110 = t93 * t95;
t82 = t96 * t110 - t94 * t98;
t73 = atan2(-t74, t82);
t70 = sin(t73);
t71 = cos(t73);
t64 = -t70 * t74 + t71 * t82;
t63 = 0.1e1 / t64 ^ 2;
t109 = t93 * t98;
t103 = t100 * t99;
t107 = t97 * t96;
t87 = -t94 * t107 + t103;
t78 = -t97 * t109 + t87 * t95;
t116 = t63 * t78;
t79 = t97 * t110 + t87 * t98;
t86 = t94 * t106 + t104;
t92 = qJ(4) + pkin(12) + qJ(6);
t90 = sin(t92);
t91 = cos(t92);
t69 = t79 * t91 + t86 * t90;
t67 = 0.1e1 / t69 ^ 2;
t68 = t79 * t90 - t86 * t91;
t115 = t67 * t68;
t114 = t71 * t74;
t81 = 0.1e1 / t82 ^ 2;
t113 = t74 * t81;
t112 = t78 ^ 2 * t63;
t111 = t86 * t98;
t108 = t93 * t99;
t102 = t68 ^ 2 * t67 + 0.1e1;
t76 = -t95 * t105 + t85 * t98;
t101 = -t70 * t82 - t114;
t84 = t94 * t103 - t107;
t83 = t96 * t109 + t94 * t95;
t80 = 0.1e1 / t82;
t72 = 0.1e1 / (t74 ^ 2 * t81 + 0.1e1);
t66 = 0.1e1 / t69;
t65 = 0.1e1 / t102;
t62 = 0.1e1 / t64;
t61 = 0.1e1 / (0.1e1 + t112);
t60 = (t108 * t113 - t80 * t84) * t95 * t72;
t59 = (t83 * t113 - t76 * t80) * t72;
t58 = t102 * t65;
t1 = [-t78 * t80 * t72, t60, t59, 0, 0, 0; (-t74 * t62 - (-t70 + (t80 * t114 + t70) * t72) * t112) * t61 (-t86 * t95 * t62 - ((t71 * t108 - t70 * t84) * t95 + t101 * t60) * t116) * t61 (t79 * t62 - (t101 * t59 - t70 * t76 + t71 * t83) * t116) * t61, 0, 0, 0; ((-t76 * t90 - t84 * t91) * t66 - (-t76 * t91 + t84 * t90) * t115) * t65 ((-t90 * t111 - t87 * t91) * t66 - (-t91 * t111 + t87 * t90) * t115) * t65 (t91 * t115 - t90 * t66) * t78 * t65, t58, 0, t58;];
Ja_rot  = t1;
