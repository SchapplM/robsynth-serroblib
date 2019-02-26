% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:45
% EndTime: 2019-02-26 19:48:45
% DurationCPUTime: 0.15s
% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
t70 = sin(pkin(10));
t72 = cos(pkin(10));
t77 = cos(qJ(2));
t73 = cos(pkin(6));
t75 = sin(qJ(2));
t81 = t73 * t75;
t63 = t70 * t77 + t72 * t81;
t69 = pkin(11) + qJ(4);
t67 = sin(t69);
t68 = cos(t69);
t71 = sin(pkin(6));
t84 = t71 * t72;
t54 = t63 * t68 - t67 * t84;
t83 = t71 * t75;
t61 = t67 * t73 + t68 * t83;
t52 = atan2(-t54, t61);
t47 = sin(t52);
t48 = cos(t52);
t43 = -t47 * t54 + t48 * t61;
t42 = 0.1e1 / t43 ^ 2;
t65 = -t70 * t81 + t72 * t77;
t85 = t70 * t71;
t57 = t65 * t68 + t67 * t85;
t90 = t42 * t57;
t56 = t65 * t67 - t68 * t85;
t74 = sin(qJ(6));
t80 = t73 * t77;
t64 = t70 * t80 + t72 * t75;
t76 = cos(qJ(6));
t86 = t64 * t76;
t50 = t56 * t74 + t86;
t46 = 0.1e1 / t50 ^ 2;
t87 = t64 * t74;
t49 = -t56 * t76 + t87;
t89 = t46 * t49;
t59 = 0.1e1 / t61 ^ 2;
t88 = t54 * t59;
t82 = t71 * t77;
t79 = t46 * t49 ^ 2 + 0.1e1;
t78 = -t47 * t61 - t48 * t54;
t62 = -t70 * t75 + t72 * t80;
t60 = -t67 * t83 + t68 * t73;
t58 = 0.1e1 / t61;
t53 = t63 * t67 + t68 * t84;
t51 = 0.1e1 / (t54 ^ 2 * t59 + 0.1e1);
t45 = 0.1e1 / t50;
t44 = 0.1e1 / t79;
t41 = 0.1e1 / t43;
t40 = 0.1e1 / (t42 * t57 ^ 2 + 0.1e1);
t39 = (-t58 * t62 + t82 * t88) * t68 * t51;
t38 = (t53 * t58 + t60 * t88) * t51;
t1 = [0, t39, 0, t38, 0, 0; 0 (-t64 * t68 * t41 - ((-t47 * t62 + t48 * t82) * t68 + t78 * t39) * t90) * t40, 0 (-t56 * t41 - (t38 * t78 + t47 * t53 + t48 * t60) * t90) * t40, 0, 0; 0 ((t65 * t74 + t67 * t86) * t45 - (t65 * t76 - t67 * t87) * t89) * t44, 0 (-t45 * t76 - t74 * t89) * t57 * t44, 0, t79 * t44;];
Ja_rot  = t1;
