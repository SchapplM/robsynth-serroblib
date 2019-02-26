% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:43
% EndTime: 2019-02-26 20:11:43
% DurationCPUTime: 0.19s
% Computational Cost: add. (948->31), mult. (1387->78), div. (95->9), fcn. (1976->13), ass. (0->55)
t95 = sin(pkin(6));
t96 = cos(pkin(11));
t110 = t95 * t96;
t101 = cos(qJ(2));
t94 = sin(pkin(11));
t105 = t94 * t101;
t97 = cos(pkin(6));
t99 = sin(qJ(2));
t108 = t97 * t99;
t87 = t96 * t108 + t105;
t93 = qJ(3) + qJ(4);
t91 = sin(t93);
t92 = cos(t93);
t78 = -t91 * t110 + t87 * t92;
t109 = t95 * t99;
t85 = t92 * t109 + t97 * t91;
t76 = atan2(-t78, t85);
t73 = sin(t76);
t74 = cos(t76);
t67 = -t73 * t78 + t74 * t85;
t66 = 0.1e1 / t67 ^ 2;
t111 = t94 * t95;
t104 = t96 * t101;
t89 = -t94 * t108 + t104;
t81 = t91 * t111 + t89 * t92;
t115 = t66 * t81;
t100 = cos(qJ(6));
t88 = t97 * t105 + t96 * t99;
t106 = t88 * t100;
t80 = -t92 * t111 + t89 * t91;
t98 = sin(qJ(6));
t72 = t80 * t98 + t106;
t70 = 0.1e1 / t72 ^ 2;
t112 = t88 * t98;
t71 = -t80 * t100 + t112;
t114 = t70 * t71;
t83 = 0.1e1 / t85 ^ 2;
t113 = t78 * t83;
t107 = t101 * t95;
t103 = t71 ^ 2 * t70 + 0.1e1;
t102 = -t73 * t85 - t74 * t78;
t86 = t97 * t104 - t94 * t99;
t84 = -t91 * t109 + t97 * t92;
t82 = 0.1e1 / t85;
t77 = t92 * t110 + t87 * t91;
t75 = 0.1e1 / (t78 ^ 2 * t83 + 0.1e1);
t69 = 0.1e1 / t72;
t68 = 0.1e1 / t103;
t65 = 0.1e1 / t67;
t64 = 0.1e1 / (t81 ^ 2 * t66 + 0.1e1);
t63 = (t107 * t113 - t82 * t86) * t92 * t75;
t62 = (t84 * t113 + t77 * t82) * t75;
t61 = (-t100 * t69 - t98 * t114) * t81 * t68;
t60 = (-t80 * t65 - (t102 * t62 + t73 * t77 + t74 * t84) * t115) * t64;
t1 = [0, t63, t62, t62, 0, 0; 0 (-t88 * t92 * t65 - ((t74 * t107 - t73 * t86) * t92 + t102 * t63) * t115) * t64, t60, t60, 0, 0; 0 ((t91 * t106 + t89 * t98) * t69 - (t89 * t100 - t91 * t112) * t114) * t68, t61, t61, 0, t103 * t68;];
Ja_rot  = t1;
