% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:55:23
% EndTime: 2019-02-26 19:55:23
% DurationCPUTime: 0.12s
% Computational Cost: add. (748->32), mult. (1032->79), div. (70->9), fcn. (1461->13), ass. (0->54)
t86 = sin(pkin(11));
t88 = cos(pkin(11));
t91 = cos(qJ(2));
t89 = cos(pkin(6));
t90 = sin(qJ(2));
t95 = t89 * t90;
t76 = t86 * t91 + t88 * t95;
t84 = pkin(12) + qJ(4);
t80 = sin(t84);
t81 = cos(t84);
t87 = sin(pkin(6));
t98 = t87 * t88;
t66 = t76 * t80 + t81 * t98;
t97 = t87 * t90;
t73 = t80 * t97 - t89 * t81;
t65 = atan2(-t66, t73);
t62 = sin(t65);
t63 = cos(t65);
t56 = -t62 * t66 + t63 * t73;
t55 = 0.1e1 / t56 ^ 2;
t78 = -t86 * t95 + t88 * t91;
t99 = t86 * t87;
t69 = t78 * t80 - t81 * t99;
t104 = t55 * t69;
t94 = t89 * t91;
t77 = t86 * t94 + t88 * t90;
t85 = qJ(5) + qJ(6);
t82 = sin(t85);
t101 = t77 * t82;
t70 = t78 * t81 + t80 * t99;
t83 = cos(t85);
t61 = t70 * t83 + t101;
t59 = 0.1e1 / t61 ^ 2;
t100 = t77 * t83;
t60 = t70 * t82 - t100;
t103 = t59 * t60;
t72 = 0.1e1 / t73 ^ 2;
t102 = t66 * t72;
t96 = t87 * t91;
t93 = t60 ^ 2 * t59 + 0.1e1;
t92 = -t62 * t73 - t63 * t66;
t75 = -t86 * t90 + t88 * t94;
t74 = t89 * t80 + t81 * t97;
t71 = 0.1e1 / t73;
t68 = t76 * t81 - t80 * t98;
t64 = 0.1e1 / (t66 ^ 2 * t72 + 0.1e1);
t58 = 0.1e1 / t61;
t57 = 0.1e1 / t93;
t54 = 0.1e1 / t56;
t53 = 0.1e1 / (t69 ^ 2 * t55 + 0.1e1);
t52 = (t96 * t102 - t71 * t75) * t80 * t64;
t51 = (t74 * t102 - t68 * t71) * t64;
t50 = t93 * t57;
t1 = [0, t52, 0, t51, 0, 0; 0 (-t77 * t80 * t54 - ((-t62 * t75 + t63 * t96) * t80 + t92 * t52) * t104) * t53, 0 (t70 * t54 - (t92 * t51 - t62 * t68 + t63 * t74) * t104) * t53, 0, 0; 0 ((-t81 * t101 - t78 * t83) * t58 - (-t81 * t100 + t78 * t82) * t103) * t57, 0 (t83 * t103 - t58 * t82) * t69 * t57, t50, t50;];
Ja_rot  = t1;
