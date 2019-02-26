% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR6_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:05
% EndTime: 2019-02-26 20:07:05
% DurationCPUTime: 0.24s
% Computational Cost: add. (555->40), mult. (1679->96), div. (60->9), fcn. (2302->15), ass. (0->59)
t94 = cos(pkin(6));
t98 = cos(qJ(2));
t106 = t94 * t98;
t88 = sin(pkin(12));
t92 = cos(pkin(12));
t96 = sin(qJ(2));
t100 = t92 * t106 - t88 * t96;
t89 = sin(pkin(7));
t90 = sin(pkin(6));
t109 = t90 * t89;
t93 = cos(pkin(7));
t116 = t100 * t93 - t92 * t109;
t107 = t94 * t96;
t82 = t92 * t107 + t88 * t98;
t95 = sin(qJ(3));
t97 = cos(qJ(3));
t67 = -t116 * t97 + t82 * t95;
t108 = t93 * t97;
t110 = t89 * t94;
t76 = -t97 * t110 + (-t98 * t108 + t95 * t96) * t90;
t66 = atan2(-t67, t76);
t63 = sin(t66);
t64 = cos(t66);
t57 = -t63 * t67 + t64 * t76;
t56 = 0.1e1 / t57 ^ 2;
t103 = t88 * t109;
t84 = -t88 * t107 + t92 * t98;
t111 = t84 * t95;
t83 = -t88 * t106 - t92 * t96;
t70 = -t97 * t103 - t83 * t108 + t111;
t115 = t56 * t70;
t71 = t84 * t97 + (t83 * t93 + t103) * t95;
t78 = t88 * t90 * t93 - t83 * t89;
t87 = sin(pkin(13));
t91 = cos(pkin(13));
t62 = t71 * t91 + t78 * t87;
t60 = 0.1e1 / t62 ^ 2;
t61 = t71 * t87 - t78 * t91;
t114 = t60 * t61;
t75 = 0.1e1 / t76 ^ 2;
t113 = t67 * t75;
t112 = t84 * t89;
t105 = t95 * t98;
t104 = t96 * t97;
t101 = -t63 * t76 - t64 * t67;
t81 = (t93 * t104 + t105) * t90;
t77 = t95 * t110 + (t93 * t105 + t104) * t90;
t74 = 0.1e1 / t76;
t73 = -t93 * t111 + t83 * t97;
t72 = t100 * t95 + t82 * t108;
t69 = t116 * t95 + t82 * t97;
t65 = 0.1e1 / (t67 ^ 2 * t75 + 0.1e1);
t59 = 0.1e1 / t62;
t58 = 0.1e1 / (t61 ^ 2 * t60 + 0.1e1);
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (t70 ^ 2 * t56 + 0.1e1);
t53 = (t81 * t113 - t72 * t74) * t65;
t52 = (t77 * t113 - t69 * t74) * t65;
t1 = [0, t53, t52, 0, 0, 0; 0 ((t84 * t108 + t83 * t95) * t55 - (t101 * t53 - t63 * t72 + t64 * t81) * t115) * t54 (t71 * t55 - (t101 * t52 - t63 * t69 + t64 * t77) * t115) * t54, 0, 0, 0; 0 ((-t112 * t91 + t73 * t87) * t59 - (t112 * t87 + t73 * t91) * t114) * t58 (t114 * t91 - t87 * t59) * t70 * t58, 0, 0, 0;];
Ja_rot  = t1;
