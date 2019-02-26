% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:46
% EndTime: 2019-02-26 19:58:46
% DurationCPUTime: 0.12s
% Computational Cost: add. (523->27), mult. (731->69), div. (57->9), fcn. (1053->11), ass. (0->44)
t65 = sin(pkin(10));
t67 = cos(pkin(10));
t70 = cos(qJ(2));
t68 = cos(pkin(6));
t69 = sin(qJ(2));
t73 = t68 * t69;
t59 = t65 * t70 + t67 * t73;
t64 = qJ(3) + pkin(11);
t62 = sin(t64);
t63 = cos(t64);
t66 = sin(pkin(6));
t76 = t66 * t67;
t47 = t59 * t62 + t63 * t76;
t75 = t66 * t69;
t54 = t62 * t75 - t68 * t63;
t45 = atan2(-t47, t54);
t42 = sin(t45);
t43 = cos(t45);
t41 = -t42 * t47 + t43 * t54;
t40 = 0.1e1 / t41 ^ 2;
t61 = -t65 * t73 + t67 * t70;
t77 = t65 * t66;
t50 = t61 * t62 - t63 * t77;
t79 = t40 * t50;
t53 = 0.1e1 / t54 ^ 2;
t78 = t47 * t53;
t74 = t66 * t70;
t72 = t68 * t70;
t71 = -t42 * t54 - t43 * t47;
t60 = t65 * t72 + t67 * t69;
t58 = -t65 * t69 + t67 * t72;
t57 = 0.1e1 / t60 ^ 2;
t56 = 0.1e1 / t60;
t55 = t68 * t62 + t63 * t75;
t52 = 0.1e1 / t54;
t51 = t61 * t63 + t62 * t77;
t49 = t59 * t63 - t62 * t76;
t46 = 0.1e1 / (t51 ^ 2 * t57 + 0.1e1);
t44 = 0.1e1 / (t47 ^ 2 * t53 + 0.1e1);
t39 = 0.1e1 / t41;
t38 = 0.1e1 / (t50 ^ 2 * t40 + 0.1e1);
t37 = (-t52 * t58 + t74 * t78) * t62 * t44;
t36 = (-t49 * t52 + t55 * t78) * t44;
t1 = [0, t37, t36, 0, 0, 0; 0 (-t60 * t62 * t39 - ((-t42 * t58 + t43 * t74) * t62 + t71 * t37) * t79) * t38 (t51 * t39 - (t71 * t36 - t42 * t49 + t43 * t55) * t79) * t38, 0, 0, 0; 0 (-t51 * t57 * t61 - t56 * t60 * t63) * t46, -t50 * t56 * t46, 0, 0, 0;];
Ja_rot  = t1;
