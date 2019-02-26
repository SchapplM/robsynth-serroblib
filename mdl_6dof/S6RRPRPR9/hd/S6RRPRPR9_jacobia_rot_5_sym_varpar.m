% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:33
% EndTime: 2019-02-26 21:42:33
% DurationCPUTime: 0.17s
% Computational Cost: add. (833->38), mult. (1217->89), div. (80->9), fcn. (1746->13), ass. (0->55)
t80 = cos(pkin(6));
t81 = sin(qJ(2));
t84 = cos(qJ(1));
t87 = t84 * t81;
t82 = sin(qJ(1));
t83 = cos(qJ(2));
t88 = t82 * t83;
t69 = t80 * t87 + t88;
t76 = pkin(11) + qJ(4);
t74 = sin(t76);
t75 = cos(t76);
t78 = sin(pkin(6));
t90 = t78 * t84;
t58 = t69 * t74 + t75 * t90;
t93 = t78 * t81;
t66 = t74 * t93 - t80 * t75;
t57 = atan2(-t58, t66);
t52 = sin(t57);
t53 = cos(t57);
t48 = -t52 * t58 + t53 * t66;
t47 = 0.1e1 / t48 ^ 2;
t86 = t84 * t83;
t89 = t82 * t81;
t71 = -t80 * t89 + t86;
t92 = t78 * t82;
t62 = t71 * t74 - t75 * t92;
t100 = t47 * t62;
t63 = t71 * t75 + t74 * t92;
t79 = cos(pkin(12));
t70 = t80 * t88 + t87;
t77 = sin(pkin(12));
t95 = t70 * t77;
t55 = t63 * t79 + t95;
t51 = 0.1e1 / t55 ^ 2;
t94 = t70 * t79;
t54 = t63 * t77 - t94;
t99 = t51 * t54;
t98 = t53 * t58;
t65 = 0.1e1 / t66 ^ 2;
t97 = t58 * t65;
t96 = t62 ^ 2 * t47;
t91 = t78 * t83;
t60 = t69 * t75 - t74 * t90;
t85 = -t52 * t66 - t98;
t68 = t80 * t86 - t89;
t67 = t80 * t74 + t75 * t93;
t64 = 0.1e1 / t66;
t56 = 0.1e1 / (t58 ^ 2 * t65 + 0.1e1);
t50 = 0.1e1 / t55;
t49 = 0.1e1 / (t54 ^ 2 * t51 + 0.1e1);
t46 = 0.1e1 / t48;
t45 = 0.1e1 / (0.1e1 + t96);
t44 = (-t64 * t68 + t91 * t97) * t74 * t56;
t43 = (-t60 * t64 + t67 * t97) * t56;
t1 = [-t62 * t64 * t56, t44, 0, t43, 0, 0; (-t58 * t46 - (-t52 + (t64 * t98 + t52) * t56) * t96) * t45 (-t70 * t74 * t46 - ((-t52 * t68 + t53 * t91) * t74 + t85 * t44) * t100) * t45, 0 (t63 * t46 - (t85 * t43 - t52 * t60 + t53 * t67) * t100) * t45, 0, 0; ((-t60 * t77 - t68 * t79) * t50 - (-t60 * t79 + t68 * t77) * t99) * t49 ((-t71 * t79 - t75 * t95) * t50 - (t71 * t77 - t75 * t94) * t99) * t49, 0 (-t77 * t50 + t79 * t99) * t62 * t49, 0, 0;];
Ja_rot  = t1;
