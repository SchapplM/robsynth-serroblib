% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:08
% EndTime: 2019-02-26 21:40:08
% DurationCPUTime: 0.26s
% Computational Cost: add. (808->39), mult. (2210->96), div. (80->9), fcn. (3119->15), ass. (0->55)
t100 = cos(qJ(4));
t102 = cos(qJ(1));
t93 = sin(pkin(6));
t108 = t102 * t93;
t101 = cos(qJ(2));
t92 = sin(pkin(11));
t95 = cos(pkin(11));
t98 = sin(qJ(2));
t106 = t101 * t95 - t98 * t92;
t105 = t101 * t92 + t98 * t95;
t96 = cos(pkin(6));
t85 = t105 * t96;
t99 = sin(qJ(1));
t74 = t102 * t85 + t106 * t99;
t97 = sin(qJ(4));
t67 = t100 * t108 + t74 * t97;
t84 = t105 * t93;
t80 = -t96 * t100 + t84 * t97;
t66 = atan2(-t67, t80);
t63 = sin(t66);
t64 = cos(t66);
t57 = -t63 * t67 + t64 * t80;
t56 = 0.1e1 / t57 ^ 2;
t104 = t102 * t106 - t99 * t85;
t110 = t93 * t99;
t71 = -t100 * t110 + t104 * t97;
t115 = t56 * t71;
t72 = t100 * t104 + t97 * t110;
t103 = t106 * t96;
t76 = -t102 * t105 - t99 * t103;
t91 = sin(pkin(12));
t94 = cos(pkin(12));
t62 = t72 * t94 - t76 * t91;
t60 = 0.1e1 / t62 ^ 2;
t61 = t72 * t91 + t76 * t94;
t114 = t60 * t61;
t113 = t64 * t67;
t79 = 0.1e1 / t80 ^ 2;
t112 = t67 * t79;
t111 = t71 ^ 2 * t56;
t109 = t100 * t76;
t69 = t74 * t100 - t97 * t108;
t107 = -t63 * t80 - t113;
t83 = t106 * t93;
t81 = t84 * t100 + t96 * t97;
t78 = 0.1e1 / t80;
t73 = t102 * t103 - t105 * t99;
t65 = 0.1e1 / (t67 ^ 2 * t79 + 0.1e1);
t59 = 0.1e1 / t62;
t58 = 0.1e1 / (t61 ^ 2 * t60 + 0.1e1);
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (0.1e1 + t111);
t53 = (t83 * t112 - t73 * t78) * t97 * t65;
t52 = (t81 * t112 - t69 * t78) * t65;
t1 = [-t71 * t78 * t65, t53, 0, t52, 0, 0; (-t67 * t55 - (-t63 + (t78 * t113 + t63) * t65) * t111) * t54 (t76 * t97 * t55 - ((-t63 * t73 + t64 * t83) * t97 + t107 * t53) * t115) * t54, 0 (t72 * t55 - (t107 * t52 - t63 * t69 + t64 * t81) * t115) * t54, 0, 0; ((-t69 * t91 - t73 * t94) * t59 - (-t69 * t94 + t73 * t91) * t114) * t58 ((-t104 * t94 + t91 * t109) * t59 - (t104 * t91 + t94 * t109) * t114) * t58, 0 (t94 * t114 - t91 * t59) * t71 * t58, 0, 0;];
Ja_rot  = t1;
