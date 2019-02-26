% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:01
% EndTime: 2019-02-26 19:52:01
% DurationCPUTime: 0.23s
% Computational Cost: add. (1149->41), mult. (1993->102), div. (87->9), fcn. (2807->13), ass. (0->60)
t100 = cos(qJ(5));
t101 = cos(qJ(2));
t96 = cos(pkin(10));
t104 = t96 * t101;
t94 = sin(pkin(10));
t97 = cos(pkin(6));
t99 = sin(qJ(2));
t102 = t97 * t104 - t94 * t99;
t95 = sin(pkin(6));
t110 = t95 * t96;
t105 = t94 * t101;
t108 = t97 * t99;
t88 = t96 * t108 + t105;
t93 = pkin(11) + qJ(4);
t91 = sin(t93);
t92 = cos(t93);
t77 = -t91 * t110 + t88 * t92;
t98 = sin(qJ(5));
t69 = t102 * t100 + t77 * t98;
t109 = t95 * t99;
t86 = t92 * t109 + t97 * t91;
t82 = t95 * t101 * t100 + t86 * t98;
t66 = atan2(-t69, t82);
t63 = sin(t66);
t64 = cos(t66);
t61 = -t63 * t69 + t64 * t82;
t60 = 0.1e1 / t61 ^ 2;
t89 = t97 * t105 + t96 * t99;
t106 = t89 * t100;
t111 = t94 * t95;
t90 = -t94 * t108 + t104;
t79 = t91 * t111 + t90 * t92;
t72 = t79 * t98 - t106;
t116 = t60 * t72;
t73 = t79 * t100 + t89 * t98;
t68 = 0.1e1 / t73 ^ 2;
t78 = t92 * t111 - t90 * t91;
t115 = t68 * t78;
t81 = 0.1e1 / t82 ^ 2;
t114 = t69 * t81;
t113 = t78 ^ 2 * t68;
t112 = t92 * t98;
t107 = t101 * t98;
t103 = -t63 * t82 - t64 * t69;
t85 = -t91 * t109 + t97 * t92;
t84 = (-t100 * t99 + t92 * t107) * t95;
t83 = t86 * t100 - t95 * t107;
t80 = 0.1e1 / t82;
t76 = -t92 * t110 - t88 * t91;
t74 = -t88 * t100 + t102 * t112;
t71 = t77 * t100 - t102 * t98;
t67 = 0.1e1 / t73;
t65 = 0.1e1 / (t69 ^ 2 * t81 + 0.1e1);
t62 = 0.1e1 / (0.1e1 + t113);
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (t72 ^ 2 * t60 + 0.1e1);
t57 = (t85 * t114 - t76 * t80) * t98 * t65;
t56 = (t84 * t114 - t74 * t80) * t65;
t55 = (t83 * t114 - t71 * t80) * t65;
t1 = [0, t56, 0, t57, t55, 0; 0 ((-t90 * t100 - t89 * t112) * t59 - (t103 * t56 - t63 * t74 + t64 * t84) * t116) * t58, 0 (t78 * t98 * t59 - ((-t63 * t76 + t64 * t85) * t98 + t103 * t57) * t116) * t58 (t73 * t59 - (t103 * t55 - t63 * t71 + t64 * t83) * t116) * t58, 0; 0 (t89 * t91 * t67 - (-t106 * t92 + t90 * t98) * t115) * t62, 0 (-t100 * t113 - t67 * t79) * t62, t72 * t62 * t115, 0;];
Ja_rot  = t1;
