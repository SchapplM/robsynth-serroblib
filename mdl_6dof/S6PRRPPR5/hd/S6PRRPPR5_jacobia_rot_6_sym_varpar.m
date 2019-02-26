% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:37
% EndTime: 2019-02-26 20:00:37
% DurationCPUTime: 0.14s
% Computational Cost: add. (382->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->50)
t71 = sin(pkin(10));
t73 = cos(pkin(10));
t78 = cos(qJ(2));
t74 = cos(pkin(6));
t76 = sin(qJ(2));
t82 = t74 * t76;
t62 = t71 * t78 + t73 * t82;
t77 = cos(qJ(3));
t72 = sin(pkin(6));
t75 = sin(qJ(3));
t85 = t72 * t75;
t55 = t62 * t77 - t73 * t85;
t84 = t72 * t77;
t66 = t74 * t75 + t76 * t84;
t53 = atan2(-t55, t66);
t50 = sin(t53);
t51 = cos(t53);
t44 = -t50 * t55 + t51 * t66;
t43 = 0.1e1 / t44 ^ 2;
t64 = -t71 * t82 + t73 * t78;
t58 = t64 * t77 + t71 * t85;
t89 = t43 * t58;
t57 = t64 * t75 - t71 * t84;
t81 = t74 * t78;
t63 = t71 * t81 + t73 * t76;
t70 = pkin(11) + qJ(6);
t68 = sin(t70);
t69 = cos(t70);
t49 = t57 * t68 + t63 * t69;
t47 = 0.1e1 / t49 ^ 2;
t48 = -t57 * t69 + t63 * t68;
t88 = t47 * t48;
t60 = 0.1e1 / t66 ^ 2;
t87 = t55 * t60;
t86 = t63 * t75;
t83 = t72 * t78;
t80 = t48 ^ 2 * t47 + 0.1e1;
t79 = -t50 * t66 - t51 * t55;
t65 = t74 * t77 - t76 * t85;
t61 = -t71 * t76 + t73 * t81;
t59 = 0.1e1 / t66;
t54 = t62 * t75 + t73 * t84;
t52 = 0.1e1 / (t55 ^ 2 * t60 + 0.1e1);
t46 = 0.1e1 / t49;
t45 = 0.1e1 / t80;
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (t58 ^ 2 * t43 + 0.1e1);
t40 = (-t59 * t61 + t83 * t87) * t77 * t52;
t39 = (t54 * t59 + t65 * t87) * t52;
t1 = [0, t40, t39, 0, 0, 0; 0 (-t63 * t77 * t42 - ((-t50 * t61 + t51 * t83) * t77 + t79 * t40) * t89) * t41 (-t57 * t42 - (t79 * t39 + t50 * t54 + t51 * t65) * t89) * t41, 0, 0, 0; 0 ((t64 * t68 + t69 * t86) * t46 - (t64 * t69 - t68 * t86) * t88) * t45 (-t46 * t69 - t68 * t88) * t58 * t45, 0, 0, t80 * t45;];
Ja_rot  = t1;
