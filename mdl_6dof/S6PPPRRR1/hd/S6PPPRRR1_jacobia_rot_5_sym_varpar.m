% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPPRRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobia_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:44
% EndTime: 2019-02-26 19:38:44
% DurationCPUTime: 0.21s
% Computational Cost: add. (968->40), mult. (2836->90), div. (35->9), fcn. (3836->19), ass. (0->56)
t95 = cos(pkin(12));
t98 = cos(pkin(6));
t111 = t95 * t98;
t88 = sin(pkin(13));
t89 = sin(pkin(12));
t94 = cos(pkin(13));
t107 = t94 * t111 - t89 * t88;
t92 = sin(pkin(6));
t112 = t95 * t92;
t91 = sin(pkin(7));
t97 = cos(pkin(7));
t105 = t107 * t97 - t91 * t112;
t106 = t88 * t111 + t89 * t94;
t87 = sin(pkin(14));
t90 = sin(pkin(8));
t93 = cos(pkin(14));
t96 = cos(pkin(8));
t120 = (t105 * t93 - t106 * t87) * t96 + (-t107 * t91 - t97 * t112) * t90;
t114 = t92 * t91;
t116 = t89 * t98;
t85 = -t94 * t116 - t95 * t88;
t108 = t89 * t114 + t85 * t97;
t86 = -t88 * t116 + t95 * t94;
t77 = t108 * t93 - t86 * t87;
t83 = t89 * t92 * t97 - t85 * t91;
t119 = t77 * t96 + t83 * t90;
t115 = t91 * t98;
t113 = t94 * t97;
t101 = cos(qJ(5));
t100 = sin(qJ(4));
t102 = cos(qJ(4));
t78 = t108 * t87 + t86 * t93;
t69 = t119 * t100 + t78 * t102;
t73 = -t77 * t90 + t83 * t96;
t99 = sin(qJ(5));
t60 = t69 * t101 + t73 * t99;
t58 = 0.1e1 / t60 ^ 2;
t59 = -t73 * t101 + t69 * t99;
t110 = t59 ^ 2 * t58 + 0.1e1;
t109 = (t93 * t115 + (t93 * t113 - t87 * t88) * t92) * t96 + (-t94 * t114 + t98 * t97) * t90;
t82 = t92 * t88 * t93 + (t92 * t113 + t115) * t87;
t76 = t105 * t87 + t106 * t93;
t72 = t109 * t100 + t82 * t102;
t71 = t82 * t100 - t109 * t102;
t70 = 0.1e1 / t71 ^ 2;
t68 = t78 * t100 - t119 * t102;
t67 = t120 * t100 + t76 * t102;
t65 = t76 * t100 - t120 * t102;
t64 = atan2(-t65, t71);
t62 = cos(t64);
t61 = sin(t64);
t57 = 0.1e1 / t110;
t56 = -t61 * t65 + t62 * t71;
t55 = 0.1e1 / t56 ^ 2;
t53 = (-t67 / t71 + t72 * t65 * t70) / (t65 ^ 2 * t70 + 0.1e1);
t1 = [0, 0, 0, t53, 0, 0; 0, 0, 0 (t69 / t56 - (-t61 * t67 + t62 * t72 + (-t61 * t71 - t62 * t65) * t53) * t68 * t55) / (t68 ^ 2 * t55 + 0.1e1) 0, 0; 0, 0, 0 (-t99 / t60 + t101 * t59 * t58) * t68 * t57, t110 * t57, 0;];
Ja_rot  = t1;
