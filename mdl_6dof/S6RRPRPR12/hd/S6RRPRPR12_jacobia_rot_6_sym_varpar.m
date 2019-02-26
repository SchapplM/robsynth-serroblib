% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:22
% EndTime: 2019-02-26 21:44:22
% DurationCPUTime: 0.16s
% Computational Cost: add. (878->37), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
t81 = cos(pkin(6));
t86 = cos(qJ(2));
t87 = cos(qJ(1));
t92 = t87 * t86;
t83 = sin(qJ(2));
t84 = sin(qJ(1));
t95 = t84 * t83;
t73 = -t81 * t92 + t95;
t79 = qJ(4) + pkin(11);
t77 = sin(t79);
t78 = cos(t79);
t80 = sin(pkin(6));
t96 = t80 * t87;
t65 = t73 * t78 + t77 * t96;
t97 = t80 * t86;
t71 = t81 * t77 + t78 * t97;
t62 = atan2(t65, t71);
t55 = sin(t62);
t56 = cos(t62);
t53 = t55 * t65 + t56 * t71;
t52 = 0.1e1 / t53 ^ 2;
t93 = t87 * t83;
t94 = t84 * t86;
t88 = t81 * t94 + t93;
t98 = t80 * t84;
t63 = t77 * t98 - t88 * t78;
t106 = t52 * t63;
t105 = t56 * t65;
t75 = -t81 * t95 + t92;
t82 = sin(qJ(6));
t101 = t75 * t82;
t64 = t88 * t77 + t78 * t98;
t85 = cos(qJ(6));
t61 = t64 * t85 + t101;
t58 = 0.1e1 / t61 ^ 2;
t100 = t75 * t85;
t60 = t64 * t82 - t100;
t104 = t58 * t60;
t103 = t63 ^ 2 * t52;
t70 = 0.1e1 / t71 ^ 2;
t102 = t65 * t70;
t99 = t80 * t83;
t91 = t60 ^ 2 * t58 + 0.1e1;
t90 = -t55 * t71 + t105;
t89 = -t73 * t77 + t78 * t96;
t74 = t81 * t93 + t94;
t72 = -t77 * t97 + t81 * t78;
t69 = 0.1e1 / t71;
t59 = 0.1e1 / (t65 ^ 2 * t70 + 0.1e1);
t57 = 0.1e1 / t61;
t54 = 0.1e1 / t91;
t51 = 0.1e1 / t53;
t50 = 0.1e1 / (0.1e1 + t103);
t49 = (t99 * t102 + t69 * t74) * t78 * t59;
t48 = (-t72 * t102 + t69 * t89) * t59;
t1 = [-t63 * t69 * t59, t49, 0, t48, 0, 0; (t65 * t51 - (-t55 + (-t69 * t105 + t55) * t59) * t103) * t50 (-t75 * t78 * t51 - ((t55 * t74 - t56 * t99) * t78 + t90 * t49) * t106) * t50, 0 (t64 * t51 - (t90 * t48 + t55 * t89 + t56 * t72) * t106) * t50, 0, 0; ((t74 * t85 + t82 * t89) * t57 - (-t74 * t82 + t85 * t89) * t104) * t54 ((t77 * t101 + t88 * t85) * t57 - (t77 * t100 - t88 * t82) * t104) * t54, 0 (t85 * t104 - t82 * t57) * t63 * t54, 0, t91 * t54;];
Ja_rot  = t1;
