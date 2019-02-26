% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:13
% EndTime: 2019-02-26 19:58:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (689->32), mult. (949->80), div. (65->9), fcn. (1349->13), ass. (0->52)
t80 = sin(pkin(10));
t82 = cos(pkin(10));
t85 = cos(qJ(2));
t83 = cos(pkin(6));
t84 = sin(qJ(2));
t89 = t83 * t84;
t70 = t80 * t85 + t82 * t89;
t79 = qJ(3) + pkin(11);
t75 = sin(t79);
t77 = cos(t79);
t81 = sin(pkin(6));
t92 = t81 * t82;
t60 = t70 * t75 + t77 * t92;
t91 = t81 * t84;
t67 = t75 * t91 - t83 * t77;
t59 = atan2(-t60, t67);
t56 = sin(t59);
t57 = cos(t59);
t50 = -t56 * t60 + t57 * t67;
t49 = 0.1e1 / t50 ^ 2;
t72 = -t80 * t89 + t82 * t85;
t93 = t80 * t81;
t63 = t72 * t75 - t77 * t93;
t97 = t49 * t63;
t64 = t72 * t77 + t75 * t93;
t88 = t83 * t85;
t71 = t80 * t88 + t82 * t84;
t78 = pkin(12) + qJ(6);
t74 = sin(t78);
t76 = cos(t78);
t55 = t64 * t76 + t71 * t74;
t53 = 0.1e1 / t55 ^ 2;
t54 = t64 * t74 - t71 * t76;
t96 = t53 * t54;
t66 = 0.1e1 / t67 ^ 2;
t95 = t60 * t66;
t94 = t71 * t77;
t90 = t81 * t85;
t87 = t54 ^ 2 * t53 + 0.1e1;
t86 = -t56 * t67 - t57 * t60;
t69 = -t80 * t84 + t82 * t88;
t68 = t83 * t75 + t77 * t91;
t65 = 0.1e1 / t67;
t62 = t70 * t77 - t75 * t92;
t58 = 0.1e1 / (t60 ^ 2 * t66 + 0.1e1);
t52 = 0.1e1 / t55;
t51 = 0.1e1 / t87;
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (t63 ^ 2 * t49 + 0.1e1);
t46 = (-t65 * t69 + t90 * t95) * t75 * t58;
t45 = (-t62 * t65 + t68 * t95) * t58;
t1 = [0, t46, t45, 0, 0, 0; 0 (-t71 * t75 * t48 - ((-t56 * t69 + t57 * t90) * t75 + t86 * t46) * t97) * t47 (t64 * t48 - (t86 * t45 - t56 * t62 + t57 * t68) * t97) * t47, 0, 0, 0; 0 ((-t72 * t76 - t74 * t94) * t52 - (t72 * t74 - t76 * t94) * t96) * t51 (-t52 * t74 + t76 * t96) * t63 * t51, 0, 0, t87 * t51;];
Ja_rot  = t1;
