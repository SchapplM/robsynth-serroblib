% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR11_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:21
% EndTime: 2019-02-26 22:36:21
% DurationCPUTime: 0.17s
% Computational Cost: add. (525->38), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->55)
t81 = cos(pkin(6));
t83 = sin(qJ(2));
t87 = cos(qJ(1));
t91 = t87 * t83;
t84 = sin(qJ(1));
t86 = cos(qJ(2));
t92 = t84 * t86;
t72 = t81 * t91 + t92;
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t80 = sin(pkin(6));
t94 = t80 * t87;
t61 = t72 * t82 + t85 * t94;
t97 = t80 * t82;
t69 = -t81 * t85 + t83 * t97;
t60 = atan2(-t61, t69);
t57 = sin(t60);
t58 = cos(t60);
t51 = -t57 * t61 + t58 * t69;
t50 = 0.1e1 / t51 ^ 2;
t90 = t87 * t86;
t93 = t84 * t83;
t74 = -t81 * t93 + t90;
t96 = t80 * t85;
t65 = t74 * t82 - t84 * t96;
t103 = t50 * t65;
t66 = t74 * t85 + t84 * t97;
t73 = t81 * t92 + t91;
t79 = qJ(4) + pkin(12);
t77 = sin(t79);
t78 = cos(t79);
t56 = t66 * t78 + t73 * t77;
t54 = 0.1e1 / t56 ^ 2;
t55 = t66 * t77 - t73 * t78;
t102 = t54 * t55;
t101 = t58 * t61;
t68 = 0.1e1 / t69 ^ 2;
t100 = t61 * t68;
t99 = t65 ^ 2 * t50;
t98 = t73 * t85;
t95 = t80 * t86;
t89 = t55 ^ 2 * t54 + 0.1e1;
t63 = t72 * t85 - t82 * t94;
t88 = -t57 * t69 - t101;
t71 = t81 * t90 - t93;
t70 = t81 * t82 + t83 * t96;
t67 = 0.1e1 / t69;
t59 = 0.1e1 / (t61 ^ 2 * t68 + 0.1e1);
t53 = 0.1e1 / t56;
t52 = 0.1e1 / t89;
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (0.1e1 + t99);
t47 = (t95 * t100 - t67 * t71) * t82 * t59;
t46 = (t70 * t100 - t63 * t67) * t59;
t1 = [-t65 * t67 * t59, t47, t46, 0, 0, 0; (-t61 * t49 - (-t57 + (t67 * t101 + t57) * t59) * t99) * t48 (-t73 * t82 * t49 - ((-t57 * t71 + t58 * t95) * t82 + t88 * t47) * t103) * t48 (t66 * t49 - (t88 * t46 - t57 * t63 + t58 * t70) * t103) * t48, 0, 0, 0; ((-t63 * t77 - t71 * t78) * t53 - (-t63 * t78 + t71 * t77) * t102) * t52 ((-t74 * t78 - t77 * t98) * t53 - (t74 * t77 - t78 * t98) * t102) * t52 (t78 * t102 - t77 * t53) * t65 * t52, t89 * t52, 0, 0;];
Ja_rot  = t1;
