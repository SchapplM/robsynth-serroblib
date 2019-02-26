% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:33
% EndTime: 2019-02-26 19:51:33
% DurationCPUTime: 0.13s
% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
t76 = sin(pkin(10));
t78 = cos(pkin(10));
t83 = cos(qJ(2));
t79 = cos(pkin(6));
t81 = sin(qJ(2));
t87 = t79 * t81;
t69 = t76 * t83 + t78 * t87;
t75 = pkin(11) + qJ(4);
t73 = sin(t75);
t74 = cos(t75);
t77 = sin(pkin(6));
t90 = t77 * t78;
t59 = t69 * t73 + t74 * t90;
t89 = t77 * t81;
t66 = t73 * t89 - t79 * t74;
t58 = atan2(-t59, t66);
t53 = sin(t58);
t54 = cos(t58);
t49 = -t53 * t59 + t54 * t66;
t48 = 0.1e1 / t49 ^ 2;
t71 = -t76 * t87 + t78 * t83;
t91 = t76 * t77;
t62 = t71 * t73 - t74 * t91;
t96 = t48 * t62;
t63 = t71 * t74 + t73 * t91;
t82 = cos(qJ(5));
t86 = t79 * t83;
t70 = t76 * t86 + t78 * t81;
t80 = sin(qJ(5));
t93 = t70 * t80;
t56 = t63 * t82 + t93;
t52 = 0.1e1 / t56 ^ 2;
t92 = t70 * t82;
t55 = t63 * t80 - t92;
t95 = t52 * t55;
t65 = 0.1e1 / t66 ^ 2;
t94 = t59 * t65;
t88 = t77 * t83;
t85 = t55 ^ 2 * t52 + 0.1e1;
t84 = -t53 * t66 - t54 * t59;
t68 = -t76 * t81 + t78 * t86;
t67 = t79 * t73 + t74 * t89;
t64 = 0.1e1 / t66;
t61 = t69 * t74 - t73 * t90;
t57 = 0.1e1 / (t59 ^ 2 * t65 + 0.1e1);
t51 = 0.1e1 / t56;
t50 = 0.1e1 / t85;
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (t62 ^ 2 * t48 + 0.1e1);
t45 = (-t64 * t68 + t88 * t94) * t73 * t57;
t44 = (-t61 * t64 + t67 * t94) * t57;
t1 = [0, t45, 0, t44, 0, 0; 0 (-t70 * t73 * t47 - ((-t53 * t68 + t54 * t88) * t73 + t84 * t45) * t96) * t46, 0 (t63 * t47 - (t84 * t44 - t53 * t61 + t54 * t67) * t96) * t46, 0, 0; 0 ((-t71 * t82 - t74 * t93) * t51 - (t71 * t80 - t74 * t92) * t95) * t50, 0 (-t51 * t80 + t82 * t95) * t62 * t50, t85 * t50, 0;];
Ja_rot  = t1;
