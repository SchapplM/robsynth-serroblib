% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:00
% EndTime: 2019-02-26 21:43:01
% DurationCPUTime: 0.17s
% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
t80 = cos(pkin(6));
t82 = sin(qJ(2));
t86 = cos(qJ(1));
t90 = t86 * t82;
t83 = sin(qJ(1));
t85 = cos(qJ(2));
t91 = t83 * t85;
t72 = t80 * t90 + t91;
t78 = pkin(11) + qJ(4);
t76 = sin(t78);
t77 = cos(t78);
t79 = sin(pkin(6));
t93 = t79 * t86;
t62 = t72 * t77 - t76 * t93;
t96 = t79 * t82;
t70 = t80 * t76 + t77 * t96;
t60 = atan2(-t62, t70);
t53 = sin(t60);
t54 = cos(t60);
t51 = -t53 * t62 + t54 * t70;
t50 = 0.1e1 / t51 ^ 2;
t89 = t86 * t85;
t92 = t83 * t82;
t74 = -t80 * t92 + t89;
t95 = t79 * t83;
t66 = t74 * t77 + t76 * t95;
t103 = t50 * t66;
t102 = t54 * t62;
t65 = t74 * t76 - t77 * t95;
t81 = sin(qJ(6));
t73 = t80 * t91 + t90;
t84 = cos(qJ(6));
t97 = t73 * t84;
t59 = t65 * t81 + t97;
t56 = 0.1e1 / t59 ^ 2;
t98 = t73 * t81;
t58 = -t65 * t84 + t98;
t101 = t56 * t58;
t68 = 0.1e1 / t70 ^ 2;
t100 = t62 * t68;
t99 = t66 ^ 2 * t50;
t94 = t79 * t85;
t88 = t58 ^ 2 * t56 + 0.1e1;
t87 = -t53 * t70 - t102;
t61 = t72 * t76 + t77 * t93;
t71 = t80 * t89 - t92;
t69 = -t76 * t96 + t80 * t77;
t67 = 0.1e1 / t70;
t57 = 0.1e1 / (t62 ^ 2 * t68 + 0.1e1);
t55 = 0.1e1 / t59;
t52 = 0.1e1 / t88;
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (0.1e1 + t99);
t47 = (t94 * t100 - t67 * t71) * t77 * t57;
t46 = (t69 * t100 + t61 * t67) * t57;
t1 = [-t66 * t67 * t57, t47, 0, t46, 0, 0; (-t62 * t49 - (-t53 + (t67 * t102 + t53) * t57) * t99) * t48 (-t73 * t77 * t49 - ((-t53 * t71 + t54 * t94) * t77 + t87 * t47) * t103) * t48, 0 (-t65 * t49 - (t87 * t46 + t53 * t61 + t54 * t69) * t103) * t48, 0, 0; ((t61 * t84 + t71 * t81) * t55 - (-t61 * t81 + t71 * t84) * t101) * t52 ((t74 * t81 + t76 * t97) * t55 - (t74 * t84 - t76 * t98) * t101) * t52, 0 (-t81 * t101 - t84 * t55) * t66 * t52, 0, t88 * t52;];
Ja_rot  = t1;
