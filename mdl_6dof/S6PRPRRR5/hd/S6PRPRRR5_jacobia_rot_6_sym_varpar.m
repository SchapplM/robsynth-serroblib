% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:15
% EndTime: 2019-02-26 19:56:15
% DurationCPUTime: 0.18s
% Computational Cost: add. (948->30), mult. (1387->79), div. (95->9), fcn. (1976->13), ass. (0->54)
t87 = sin(pkin(6));
t88 = cos(pkin(11));
t100 = t87 * t88;
t86 = sin(pkin(11));
t91 = sin(qJ(2));
t89 = cos(pkin(6));
t93 = cos(qJ(2));
t96 = t89 * t93;
t79 = t86 * t91 - t88 * t96;
t85 = qJ(4) + qJ(5);
t83 = sin(t85);
t84 = cos(t85);
t72 = t83 * t100 + t79 * t84;
t98 = t87 * t93;
t77 = t89 * t83 + t84 * t98;
t69 = atan2(t72, t77);
t66 = sin(t69);
t67 = cos(t69);
t60 = t66 * t72 + t67 * t77;
t59 = 0.1e1 / t60 ^ 2;
t101 = t86 * t87;
t81 = t86 * t96 + t88 * t91;
t70 = t83 * t101 - t81 * t84;
t106 = t59 * t70;
t97 = t89 * t91;
t82 = -t86 * t97 + t88 * t93;
t90 = sin(qJ(6));
t103 = t82 * t90;
t71 = t84 * t101 + t81 * t83;
t92 = cos(qJ(6));
t65 = t71 * t92 + t103;
t63 = 0.1e1 / t65 ^ 2;
t102 = t82 * t92;
t64 = t71 * t90 - t102;
t105 = t63 * t64;
t76 = 0.1e1 / t77 ^ 2;
t104 = t72 * t76;
t99 = t87 * t91;
t95 = t64 ^ 2 * t63 + 0.1e1;
t94 = -t66 * t77 + t67 * t72;
t80 = t86 * t93 + t88 * t97;
t78 = -t83 * t98 + t89 * t84;
t75 = 0.1e1 / t77;
t73 = t84 * t100 - t79 * t83;
t68 = 0.1e1 / (t72 ^ 2 * t76 + 0.1e1);
t62 = 0.1e1 / t65;
t61 = 0.1e1 / t95;
t58 = 0.1e1 / t60;
t57 = 0.1e1 / (t70 ^ 2 * t59 + 0.1e1);
t56 = (t99 * t104 + t75 * t80) * t84 * t68;
t55 = (-t78 * t104 + t73 * t75) * t68;
t54 = (t92 * t105 - t62 * t90) * t70 * t61;
t53 = (t71 * t58 - (t94 * t55 + t66 * t73 + t67 * t78) * t106) * t57;
t1 = [0, t56, 0, t55, t55, 0; 0 (-t82 * t84 * t58 - ((t66 * t80 - t67 * t99) * t84 + t94 * t56) * t106) * t57, 0, t53, t53, 0; 0 ((t83 * t103 + t81 * t92) * t62 - (t83 * t102 - t81 * t90) * t105) * t61, 0, t54, t54, t95 * t61;];
Ja_rot  = t1;
