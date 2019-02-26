% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR12_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:25
% EndTime: 2019-02-26 22:22:26
% DurationCPUTime: 0.24s
% Computational Cost: add. (525->38), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->55)
t80 = cos(pkin(6));
t82 = sin(qJ(2));
t86 = cos(qJ(1));
t90 = t86 * t82;
t83 = sin(qJ(1));
t85 = cos(qJ(2));
t91 = t83 * t85;
t71 = t80 * t90 + t91;
t81 = sin(qJ(3));
t84 = cos(qJ(3));
t79 = sin(pkin(6));
t93 = t79 * t86;
t60 = t71 * t81 + t84 * t93;
t96 = t79 * t81;
t68 = -t80 * t84 + t82 * t96;
t59 = atan2(-t60, t68);
t56 = sin(t59);
t57 = cos(t59);
t50 = -t56 * t60 + t57 * t68;
t49 = 0.1e1 / t50 ^ 2;
t89 = t86 * t85;
t92 = t83 * t82;
t73 = -t80 * t92 + t89;
t95 = t79 * t84;
t64 = t73 * t81 - t83 * t95;
t102 = t49 * t64;
t65 = t73 * t84 + t83 * t96;
t72 = t80 * t91 + t90;
t78 = pkin(12) + qJ(5);
t76 = sin(t78);
t77 = cos(t78);
t55 = t65 * t77 + t72 * t76;
t53 = 0.1e1 / t55 ^ 2;
t54 = t65 * t76 - t72 * t77;
t101 = t53 * t54;
t100 = t57 * t60;
t67 = 0.1e1 / t68 ^ 2;
t99 = t60 * t67;
t98 = t64 ^ 2 * t49;
t97 = t72 * t84;
t94 = t79 * t85;
t88 = t54 ^ 2 * t53 + 0.1e1;
t62 = t71 * t84 - t81 * t93;
t87 = -t56 * t68 - t100;
t70 = t80 * t89 - t92;
t69 = t80 * t81 + t82 * t95;
t66 = 0.1e1 / t68;
t58 = 0.1e1 / (t60 ^ 2 * t67 + 0.1e1);
t52 = 0.1e1 / t55;
t51 = 0.1e1 / t88;
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (0.1e1 + t98);
t46 = (-t66 * t70 + t94 * t99) * t81 * t58;
t45 = (-t62 * t66 + t69 * t99) * t58;
t1 = [-t64 * t66 * t58, t46, t45, 0, 0, 0; (-t60 * t48 - (-t56 + (t66 * t100 + t56) * t58) * t98) * t47 (-t72 * t81 * t48 - ((-t56 * t70 + t57 * t94) * t81 + t87 * t46) * t102) * t47 (t65 * t48 - (t87 * t45 - t56 * t62 + t57 * t69) * t102) * t47, 0, 0, 0; ((-t62 * t76 - t70 * t77) * t52 - (-t62 * t77 + t70 * t76) * t101) * t51 ((-t73 * t77 - t76 * t97) * t52 - (t73 * t76 - t77 * t97) * t101) * t51 (t77 * t101 - t76 * t52) * t64 * t51, 0, t88 * t51, 0;];
Ja_rot  = t1;
