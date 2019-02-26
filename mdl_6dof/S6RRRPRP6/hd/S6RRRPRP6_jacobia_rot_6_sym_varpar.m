% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:03
% EndTime: 2019-02-26 22:12:03
% DurationCPUTime: 0.20s
% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
t86 = sin(pkin(6));
t93 = cos(qJ(1));
t100 = t86 * t93;
t87 = cos(pkin(6));
t89 = sin(qJ(2));
t97 = t93 * t89;
t90 = sin(qJ(1));
t92 = cos(qJ(2));
t98 = t90 * t92;
t78 = t87 * t97 + t98;
t85 = qJ(3) + pkin(11);
t83 = sin(t85);
t84 = cos(t85);
t67 = t84 * t100 + t78 * t83;
t103 = t86 * t89;
t75 = t83 * t103 - t87 * t84;
t66 = atan2(-t67, t75);
t59 = sin(t66);
t60 = cos(t66);
t57 = -t59 * t67 + t60 * t75;
t56 = 0.1e1 / t57 ^ 2;
t102 = t86 * t90;
t96 = t93 * t92;
t99 = t90 * t89;
t80 = -t87 * t99 + t96;
t71 = -t84 * t102 + t80 * t83;
t110 = t56 * t71;
t109 = t60 * t67;
t79 = t87 * t98 + t97;
t88 = sin(qJ(5));
t105 = t79 * t88;
t72 = t83 * t102 + t80 * t84;
t91 = cos(qJ(5));
t65 = t72 * t91 + t105;
t62 = 0.1e1 / t65 ^ 2;
t104 = t79 * t91;
t64 = t72 * t88 - t104;
t108 = t62 * t64;
t74 = 0.1e1 / t75 ^ 2;
t107 = t67 * t74;
t106 = t71 ^ 2 * t56;
t101 = t86 * t92;
t95 = t64 ^ 2 * t62 + 0.1e1;
t69 = -t83 * t100 + t78 * t84;
t94 = -t59 * t75 - t109;
t77 = t87 * t96 - t99;
t76 = t84 * t103 + t87 * t83;
t73 = 0.1e1 / t75;
t63 = 0.1e1 / (t67 ^ 2 * t74 + 0.1e1);
t61 = 0.1e1 / t65;
t58 = 0.1e1 / t95;
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (0.1e1 + t106);
t53 = (t101 * t107 - t73 * t77) * t83 * t63;
t52 = (t76 * t107 - t69 * t73) * t63;
t1 = [-t71 * t73 * t63, t53, t52, 0, 0, 0; (-t67 * t55 - (-t59 + (t73 * t109 + t59) * t63) * t106) * t54 (-t79 * t83 * t55 - ((t60 * t101 - t59 * t77) * t83 + t94 * t53) * t110) * t54 (t72 * t55 - (t94 * t52 - t59 * t69 + t60 * t76) * t110) * t54, 0, 0, 0; ((-t69 * t88 - t77 * t91) * t61 - (-t69 * t91 + t77 * t88) * t108) * t58 ((-t84 * t105 - t80 * t91) * t61 - (-t84 * t104 + t80 * t88) * t108) * t58 (t91 * t108 - t88 * t61) * t71 * t58, 0, t95 * t58, 0;];
Ja_rot  = t1;
