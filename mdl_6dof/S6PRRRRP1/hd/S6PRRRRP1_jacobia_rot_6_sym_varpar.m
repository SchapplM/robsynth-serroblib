% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:13
% EndTime: 2019-02-26 20:15:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (948->31), mult. (1387->79), div. (95->9), fcn. (1976->13), ass. (0->54)
t101 = sin(pkin(6));
t102 = cos(pkin(11));
t114 = t101 * t102;
t100 = sin(pkin(11));
t107 = cos(qJ(2));
t103 = cos(pkin(6));
t105 = sin(qJ(2));
t111 = t103 * t105;
t93 = t100 * t107 + t102 * t111;
t99 = qJ(3) + qJ(4);
t97 = sin(t99);
t98 = cos(t99);
t83 = t98 * t114 + t93 * t97;
t113 = t101 * t105;
t90 = -t103 * t98 + t97 * t113;
t82 = atan2(-t83, t90);
t79 = sin(t82);
t80 = cos(t82);
t73 = -t79 * t83 + t80 * t90;
t72 = 0.1e1 / t73 ^ 2;
t115 = t100 * t101;
t95 = -t100 * t111 + t102 * t107;
t86 = -t98 * t115 + t95 * t97;
t120 = t72 * t86;
t106 = cos(qJ(5));
t104 = sin(qJ(5));
t110 = t103 * t107;
t94 = t100 * t110 + t102 * t105;
t117 = t94 * t104;
t87 = t97 * t115 + t95 * t98;
t78 = t87 * t106 + t117;
t76 = 0.1e1 / t78 ^ 2;
t116 = t94 * t106;
t77 = t87 * t104 - t116;
t119 = t76 * t77;
t89 = 0.1e1 / t90 ^ 2;
t118 = t83 * t89;
t112 = t101 * t107;
t109 = t77 ^ 2 * t76 + 0.1e1;
t108 = -t79 * t90 - t80 * t83;
t92 = -t100 * t105 + t102 * t110;
t91 = t103 * t97 + t98 * t113;
t88 = 0.1e1 / t90;
t85 = -t97 * t114 + t93 * t98;
t81 = 0.1e1 / (t83 ^ 2 * t89 + 0.1e1);
t75 = 0.1e1 / t78;
t74 = 0.1e1 / t109;
t71 = 0.1e1 / t73;
t70 = 0.1e1 / (t86 ^ 2 * t72 + 0.1e1);
t69 = (t112 * t118 - t88 * t92) * t97 * t81;
t68 = (t91 * t118 - t85 * t88) * t81;
t67 = (-t104 * t75 + t106 * t119) * t86 * t74;
t66 = (t87 * t71 - (t108 * t68 - t79 * t85 + t80 * t91) * t120) * t70;
t1 = [0, t69, t68, t68, 0, 0; 0 (-t94 * t97 * t71 - ((t80 * t112 - t79 * t92) * t97 + t108 * t69) * t120) * t70, t66, t66, 0, 0; 0 ((-t95 * t106 - t98 * t117) * t75 - (t95 * t104 - t98 * t116) * t119) * t74, t67, t67, t109 * t74, 0;];
Ja_rot  = t1;
