% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:03
% EndTime: 2019-02-26 19:47:03
% DurationCPUTime: 0.18s
% Computational Cost: add. (573->32), mult. (1574->83), div. (60->9), fcn. (2215->15), ass. (0->53)
t88 = cos(pkin(6));
t82 = sin(pkin(11));
t86 = cos(pkin(11));
t90 = sin(qJ(2));
t92 = cos(qJ(2));
t94 = t92 * t82 + t90 * t86;
t76 = t94 * t88;
t77 = t90 * t82 - t92 * t86;
t83 = sin(pkin(10));
t87 = cos(pkin(10));
t65 = t87 * t76 - t83 * t77;
t89 = sin(qJ(4));
t84 = sin(pkin(6));
t91 = cos(qJ(4));
t97 = t84 * t91;
t59 = t65 * t89 + t87 * t97;
t75 = t94 * t84;
t71 = t75 * t89 - t88 * t91;
t58 = atan2(-t59, t71);
t55 = sin(t58);
t56 = cos(t58);
t49 = -t55 * t59 + t56 * t71;
t48 = 0.1e1 / t49 ^ 2;
t95 = -t83 * t76 - t87 * t77;
t62 = -t83 * t97 + t89 * t95;
t102 = t48 * t62;
t98 = t84 * t89;
t63 = t83 * t98 + t91 * t95;
t93 = t77 * t88;
t67 = t83 * t93 - t87 * t94;
t81 = sin(pkin(12));
t85 = cos(pkin(12));
t54 = t63 * t85 - t67 * t81;
t52 = 0.1e1 / t54 ^ 2;
t53 = t63 * t81 + t67 * t85;
t101 = t52 * t53;
t70 = 0.1e1 / t71 ^ 2;
t100 = t59 * t70;
t99 = t67 * t91;
t96 = -t55 * t71 - t56 * t59;
t74 = t77 * t84;
t72 = t75 * t91 + t88 * t89;
t69 = 0.1e1 / t71;
t64 = -t83 * t94 - t87 * t93;
t61 = t65 * t91 - t87 * t98;
t57 = 0.1e1 / (t59 ^ 2 * t70 + 0.1e1);
t51 = 0.1e1 / t54;
t50 = 0.1e1 / (t53 ^ 2 * t52 + 0.1e1);
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (t62 ^ 2 * t48 + 0.1e1);
t45 = (-t74 * t100 - t64 * t69) * t89 * t57;
t44 = (t72 * t100 - t61 * t69) * t57;
t1 = [0, t45, 0, t44, 0, 0; 0 (t67 * t89 * t47 - ((-t55 * t64 - t56 * t74) * t89 + t96 * t45) * t102) * t46, 0 (t63 * t47 - (t96 * t44 - t55 * t61 + t56 * t72) * t102) * t46, 0, 0; 0 ((t81 * t99 - t85 * t95) * t51 - (t81 * t95 + t85 * t99) * t101) * t50, 0 (t85 * t101 - t51 * t81) * t62 * t50, 0, 0;];
Ja_rot  = t1;
