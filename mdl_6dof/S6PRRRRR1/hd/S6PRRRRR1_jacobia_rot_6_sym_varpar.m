% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR1
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
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:36
% EndTime: 2019-02-26 20:18:37
% DurationCPUTime: 0.18s
% Computational Cost: add. (1866->31), mult. (1825->77), div. (125->9), fcn. (2603->13), ass. (0->56)
t97 = sin(pkin(6));
t98 = cos(pkin(12));
t114 = t97 * t98;
t101 = sin(qJ(2));
t107 = t98 * t101;
t103 = cos(qJ(2));
t96 = sin(pkin(12));
t108 = t96 * t103;
t99 = cos(pkin(6));
t89 = t99 * t107 + t108;
t95 = qJ(3) + qJ(4) + qJ(5);
t93 = sin(t95);
t94 = cos(t95);
t79 = t94 * t114 + t89 * t93;
t113 = t101 * t97;
t86 = t93 * t113 - t99 * t94;
t76 = atan2(-t79, t86);
t71 = sin(t76);
t72 = cos(t76);
t69 = -t71 * t79 + t72 * t86;
t68 = 0.1e1 / t69 ^ 2;
t115 = t96 * t97;
t106 = t98 * t103;
t109 = t96 * t101;
t91 = -t99 * t109 + t106;
t82 = -t94 * t115 + t91 * t93;
t118 = t68 * t82;
t102 = cos(qJ(6));
t100 = sin(qJ(6));
t90 = t99 * t108 + t107;
t111 = t90 * t100;
t83 = t93 * t115 + t91 * t94;
t78 = t83 * t102 + t111;
t75 = 0.1e1 / t78 ^ 2;
t110 = t90 * t102;
t77 = t83 * t100 - t110;
t117 = t75 * t77;
t85 = 0.1e1 / t86 ^ 2;
t116 = t79 * t85;
t112 = t103 * t97;
t105 = t77 ^ 2 * t75 + 0.1e1;
t104 = -t71 * t86 - t72 * t79;
t88 = t99 * t106 - t109;
t87 = t94 * t113 + t99 * t93;
t84 = 0.1e1 / t86;
t81 = -t93 * t114 + t89 * t94;
t74 = 0.1e1 / t78;
t73 = 0.1e1 / (t79 ^ 2 * t85 + 0.1e1);
t70 = 0.1e1 / t105;
t67 = 0.1e1 / t69;
t66 = 0.1e1 / (t82 ^ 2 * t68 + 0.1e1);
t65 = (t112 * t116 - t84 * t88) * t93 * t73;
t64 = (t87 * t116 - t81 * t84) * t73;
t63 = (-t100 * t74 + t102 * t117) * t82 * t70;
t62 = (t83 * t67 - (t104 * t64 - t71 * t81 + t72 * t87) * t118) * t66;
t1 = [0, t65, t64, t64, t64, 0; 0 (-t90 * t93 * t67 - ((t72 * t112 - t71 * t88) * t93 + t104 * t65) * t118) * t66, t62, t62, t62, 0; 0 ((-t91 * t102 - t94 * t111) * t74 - (t91 * t100 - t94 * t110) * t117) * t70, t63, t63, t63, t105 * t70;];
Ja_rot  = t1;
