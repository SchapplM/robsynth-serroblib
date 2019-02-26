% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6PRRPPR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:08
% EndTime: 2019-02-26 19:58:08
% DurationCPUTime: 0.14s
% Computational Cost: add. (596->31), mult. (866->78), div. (60->9), fcn. (1237->13), ass. (0->51)
t73 = sin(pkin(10));
t76 = cos(pkin(10));
t79 = cos(qJ(2));
t77 = cos(pkin(6));
t78 = sin(qJ(2));
t82 = t77 * t78;
t65 = t73 * t79 + t76 * t82;
t71 = qJ(3) + pkin(11);
t69 = sin(t71);
t70 = cos(t71);
t74 = sin(pkin(6));
t85 = t74 * t76;
t55 = t65 * t69 + t70 * t85;
t84 = t74 * t78;
t62 = t69 * t84 - t70 * t77;
t54 = atan2(-t55, t62);
t51 = sin(t54);
t52 = cos(t54);
t45 = -t51 * t55 + t52 * t62;
t44 = 0.1e1 / t45 ^ 2;
t67 = -t73 * t82 + t76 * t79;
t86 = t73 * t74;
t58 = t67 * t69 - t70 * t86;
t91 = t44 * t58;
t59 = t67 * t70 + t69 * t86;
t75 = cos(pkin(12));
t81 = t77 * t79;
t66 = t73 * t81 + t76 * t78;
t72 = sin(pkin(12));
t88 = t66 * t72;
t50 = t59 * t75 + t88;
t48 = 0.1e1 / t50 ^ 2;
t87 = t66 * t75;
t49 = t59 * t72 - t87;
t90 = t48 * t49;
t61 = 0.1e1 / t62 ^ 2;
t89 = t55 * t61;
t83 = t74 * t79;
t80 = -t51 * t62 - t52 * t55;
t64 = -t73 * t78 + t76 * t81;
t63 = t69 * t77 + t70 * t84;
t60 = 0.1e1 / t62;
t57 = t65 * t70 - t69 * t85;
t53 = 0.1e1 / (t55 ^ 2 * t61 + 0.1e1);
t47 = 0.1e1 / t50;
t46 = 0.1e1 / (t48 * t49 ^ 2 + 0.1e1);
t43 = 0.1e1 / t45;
t42 = 0.1e1 / (t44 * t58 ^ 2 + 0.1e1);
t41 = (-t60 * t64 + t83 * t89) * t69 * t53;
t40 = (-t57 * t60 + t63 * t89) * t53;
t1 = [0, t41, t40, 0, 0, 0; 0 (-t66 * t69 * t43 - ((-t51 * t64 + t52 * t83) * t69 + t80 * t41) * t91) * t42 (t59 * t43 - (t40 * t80 - t51 * t57 + t52 * t63) * t91) * t42, 0, 0, 0; 0 ((-t67 * t75 - t70 * t88) * t47 - (t67 * t72 - t70 * t87) * t90) * t46 (-t47 * t72 + t75 * t90) * t58 * t46, 0, 0, 0;];
Ja_rot  = t1;
