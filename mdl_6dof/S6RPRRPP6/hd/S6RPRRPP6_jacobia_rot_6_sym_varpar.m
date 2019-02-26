% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:25
% EndTime: 2019-02-26 20:59:25
% DurationCPUTime: 0.13s
% Computational Cost: add. (447->25), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->41)
t59 = cos(qJ(3));
t73 = t59 ^ 2;
t57 = sin(qJ(3));
t52 = qJ(4) + pkin(9);
t50 = sin(t52);
t60 = cos(qJ(1));
t63 = t60 * t50;
t51 = cos(t52);
t58 = sin(qJ(1));
t66 = t58 * t51;
t45 = t57 * t63 + t66;
t64 = t59 * t50;
t39 = atan2(t45, t64);
t36 = sin(t39);
t37 = cos(t39);
t35 = t36 * t45 + t37 * t64;
t34 = 0.1e1 / t35 ^ 2;
t62 = t60 * t51;
t67 = t58 * t50;
t43 = t57 * t67 - t62;
t72 = t34 * t43;
t70 = t37 * t45;
t69 = t43 ^ 2 * t34;
t48 = 0.1e1 / t50;
t55 = 0.1e1 / t59;
t68 = t48 * t55;
t65 = t58 * t59;
t44 = t57 * t66 + t63;
t42 = 0.1e1 / t44 ^ 2;
t61 = t58 ^ 2 * t73 * t42;
t56 = 0.1e1 / t73;
t49 = 0.1e1 / t50 ^ 2;
t46 = t57 * t62 - t67;
t41 = 0.1e1 / t44;
t40 = 0.1e1 / (0.1e1 + t61);
t38 = 0.1e1 / (t45 ^ 2 * t56 * t49 + 0.1e1);
t33 = 0.1e1 / t35;
t32 = (t45 * t48 * t56 * t57 + t60) * t38;
t31 = 0.1e1 / (0.1e1 + t69);
t30 = (-t45 * t49 * t51 + t46 * t48) * t55 * t38;
t1 = [-t43 * t38 * t68, 0, t32, t30, 0, 0; (t45 * t33 - (-t36 + (-t68 * t70 + t36) * t38) * t69) * t31, 0 (-t32 * t70 * t72 + (t33 * t65 - (-t37 * t57 + (-t32 + t60) * t59 * t36) * t72) * t50) * t31 (t44 * t33 - (t37 * t59 * t51 + t36 * t46 + (-t36 * t64 + t70) * t30) * t72) * t31, 0, 0; (-t42 * t46 * t58 + t41 * t60) * t59 * t40, 0 (-t41 * t57 * t58 - t51 * t61) * t40, t43 * t42 * t40 * t65, 0, 0;];
Ja_rot  = t1;
