% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRP6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:15
% EndTime: 2019-02-26 19:53:15
% DurationCPUTime: 0.16s
% Computational Cost: add. (714->40), mult. (1993->102), div. (87->9), fcn. (2807->13), ass. (0->59)
t92 = cos(pkin(6));
t95 = sin(qJ(2));
t104 = t92 * t95;
t89 = sin(pkin(10));
t91 = cos(pkin(10));
t98 = cos(qJ(2));
t100 = t91 * t104 + t89 * t98;
t90 = sin(pkin(6));
t97 = cos(qJ(4));
t106 = t90 * t97;
t103 = t92 * t98;
t84 = -t91 * t103 + t89 * t95;
t94 = sin(qJ(4));
t77 = -t91 * t106 + t84 * t94;
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t69 = -t100 * t96 + t77 * t93;
t105 = t90 * t98;
t88 = -t94 * t105 + t92 * t97;
t80 = -t90 * t95 * t96 + t88 * t93;
t64 = atan2(-t69, t80);
t61 = sin(t64);
t62 = cos(t64);
t59 = -t61 * t69 + t62 * t80;
t58 = 0.1e1 / t59 ^ 2;
t86 = -t89 * t104 + t91 * t98;
t108 = t86 * t96;
t85 = t89 * t103 + t91 * t95;
t75 = t89 * t106 + t85 * t94;
t67 = t75 * t93 - t108;
t113 = t58 * t67;
t109 = t86 * t93;
t68 = t75 * t96 + t109;
t66 = 0.1e1 / t68 ^ 2;
t107 = t90 * t94;
t74 = -t89 * t107 + t85 * t97;
t112 = t66 * t74;
t79 = 0.1e1 / t80 ^ 2;
t111 = t69 * t79;
t110 = t74 ^ 2 * t66;
t102 = t93 * t95;
t101 = -t61 * t80 - t62 * t69;
t99 = t100 * t93;
t87 = -t97 * t105 - t92 * t94;
t82 = (t94 * t102 - t96 * t98) * t90;
t81 = t90 * t102 + t88 * t96;
t78 = 0.1e1 / t80;
t76 = t91 * t107 + t84 * t97;
t72 = t84 * t96 + t94 * t99;
t71 = t77 * t96 + t99;
t65 = 0.1e1 / t68;
t63 = 0.1e1 / (t69 ^ 2 * t79 + 0.1e1);
t60 = 0.1e1 / (0.1e1 + t110);
t57 = 0.1e1 / t59;
t56 = 0.1e1 / (t67 ^ 2 * t58 + 0.1e1);
t55 = (t87 * t111 - t76 * t78) * t93 * t63;
t54 = (t82 * t111 - t72 * t78) * t63;
t53 = (t81 * t111 - t71 * t78) * t63;
t1 = [0, t54, 0, t55, t53, 0; 0 ((t94 * t109 + t85 * t96) * t57 - (t101 * t54 - t61 * t72 + t62 * t82) * t113) * t56, 0 (t74 * t93 * t57 - ((-t61 * t76 + t62 * t87) * t93 + t101 * t55) * t113) * t56 (t68 * t57 - (t101 * t53 - t61 * t71 + t62 * t81) * t113) * t56, 0; 0 (t86 * t97 * t65 - (t94 * t108 - t85 * t93) * t112) * t60, 0 (-t96 * t110 - t65 * t75) * t60, t67 * t60 * t112, 0;];
Ja_rot  = t1;
