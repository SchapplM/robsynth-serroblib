% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6PRPRPR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:45
% EndTime: 2019-02-26 19:48:45
% DurationCPUTime: 0.10s
% Computational Cost: add. (523->27), mult. (731->69), div. (57->9), fcn. (1053->11), ass. (0->44)
t62 = sin(pkin(10));
t64 = cos(pkin(10));
t67 = cos(qJ(2));
t65 = cos(pkin(6));
t66 = sin(qJ(2));
t70 = t65 * t66;
t56 = t62 * t67 + t64 * t70;
t61 = pkin(11) + qJ(4);
t59 = sin(t61);
t60 = cos(t61);
t63 = sin(pkin(6));
t73 = t63 * t64;
t44 = t56 * t59 + t60 * t73;
t72 = t63 * t66;
t51 = t59 * t72 - t65 * t60;
t42 = atan2(-t44, t51);
t39 = sin(t42);
t40 = cos(t42);
t38 = -t39 * t44 + t40 * t51;
t37 = 0.1e1 / t38 ^ 2;
t58 = -t62 * t70 + t64 * t67;
t74 = t62 * t63;
t47 = t58 * t59 - t60 * t74;
t76 = t37 * t47;
t50 = 0.1e1 / t51 ^ 2;
t75 = t44 * t50;
t71 = t63 * t67;
t69 = t65 * t67;
t68 = -t39 * t51 - t40 * t44;
t57 = t62 * t69 + t64 * t66;
t55 = -t62 * t66 + t64 * t69;
t54 = 0.1e1 / t57 ^ 2;
t53 = 0.1e1 / t57;
t52 = t65 * t59 + t60 * t72;
t49 = 0.1e1 / t51;
t48 = t58 * t60 + t59 * t74;
t46 = t56 * t60 - t59 * t73;
t43 = 0.1e1 / (t48 ^ 2 * t54 + 0.1e1);
t41 = 0.1e1 / (t44 ^ 2 * t50 + 0.1e1);
t36 = 0.1e1 / t38;
t35 = 0.1e1 / (t47 ^ 2 * t37 + 0.1e1);
t34 = (-t49 * t55 + t71 * t75) * t59 * t41;
t33 = (-t46 * t49 + t52 * t75) * t41;
t1 = [0, t34, 0, t33, 0, 0; 0 (-t57 * t59 * t36 - ((-t39 * t55 + t40 * t71) * t59 + t68 * t34) * t76) * t35, 0 (t48 * t36 - (t68 * t33 - t39 * t46 + t40 * t52) * t76) * t35, 0, 0; 0 (-t48 * t54 * t58 - t53 * t57 * t60) * t43, 0, -t47 * t53 * t43, 0, 0;];
Ja_rot  = t1;
