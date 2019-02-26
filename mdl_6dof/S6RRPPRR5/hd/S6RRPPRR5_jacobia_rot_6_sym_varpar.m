% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:48
% EndTime: 2019-02-26 21:30:48
% DurationCPUTime: 0.18s
% Computational Cost: add. (459->37), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->54)
t75 = cos(pkin(6));
t78 = sin(qJ(2));
t83 = cos(qJ(1));
t87 = t83 * t78;
t79 = sin(qJ(1));
t82 = cos(qJ(2));
t88 = t79 * t82;
t68 = t75 * t87 + t88;
t77 = sin(qJ(5));
t81 = cos(qJ(5));
t74 = sin(pkin(6));
t90 = t74 * t83;
t57 = t68 * t77 - t81 * t90;
t93 = t74 * t77;
t65 = t75 * t81 + t78 * t93;
t56 = atan2(-t57, t65);
t53 = sin(t56);
t54 = cos(t56);
t47 = -t53 * t57 + t54 * t65;
t46 = 0.1e1 / t47 ^ 2;
t86 = t83 * t82;
t89 = t79 * t78;
t70 = -t75 * t89 + t86;
t92 = t74 * t81;
t61 = t70 * t77 + t79 * t92;
t99 = t46 * t61;
t62 = t70 * t81 - t79 * t93;
t69 = -t75 * t88 - t87;
t76 = sin(qJ(6));
t80 = cos(qJ(6));
t52 = t62 * t80 + t69 * t76;
t50 = 0.1e1 / t52 ^ 2;
t51 = t62 * t76 - t69 * t80;
t98 = t50 * t51;
t97 = t54 * t57;
t64 = 0.1e1 / t65 ^ 2;
t96 = t57 * t64;
t95 = t61 ^ 2 * t46;
t94 = t69 * t81;
t91 = t74 * t82;
t85 = t51 ^ 2 * t50 + 0.1e1;
t84 = -t53 * t65 - t97;
t59 = t68 * t81 + t77 * t90;
t67 = -t75 * t86 + t89;
t66 = -t75 * t77 + t78 * t92;
t63 = 0.1e1 / t65;
t55 = 0.1e1 / (t57 ^ 2 * t64 + 0.1e1);
t49 = 0.1e1 / t52;
t48 = 0.1e1 / t85;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (0.1e1 + t95);
t43 = (t63 * t67 + t91 * t96) * t77 * t55;
t42 = (-t59 * t63 + t66 * t96) * t55;
t1 = [-t61 * t63 * t55, t43, 0, 0, t42, 0; (-t57 * t45 - (-t53 + (t63 * t97 + t53) * t55) * t95) * t44 (t69 * t77 * t45 - ((t53 * t67 + t54 * t91) * t77 + t84 * t43) * t99) * t44, 0, 0 (t62 * t45 - (t84 * t42 - t53 * t59 + t54 * t66) * t99) * t44, 0; ((-t59 * t76 - t67 * t80) * t49 - (-t59 * t80 + t67 * t76) * t98) * t48 ((t70 * t80 + t76 * t94) * t49 - (-t70 * t76 + t80 * t94) * t98) * t48, 0, 0 (-t76 * t49 + t80 * t98) * t61 * t48, t85 * t48;];
Ja_rot  = t1;
