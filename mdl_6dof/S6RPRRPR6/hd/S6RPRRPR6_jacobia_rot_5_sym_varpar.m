% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:47
% EndTime: 2019-02-26 21:03:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t52 = pkin(10) + qJ(3);
t50 = cos(t52);
t48 = sin(t52);
t54 = sin(qJ(1));
t61 = t54 * t48;
t43 = atan2(-t61, -t50);
t41 = sin(t43);
t42 = cos(t43);
t34 = -t41 * t61 - t42 * t50;
t33 = 0.1e1 / t34 ^ 2;
t55 = cos(qJ(1));
t67 = t33 * t55 ^ 2;
t53 = qJ(4) + pkin(11);
t51 = cos(t53);
t57 = t55 * t51;
t49 = sin(t53);
t60 = t54 * t49;
t40 = t50 * t57 + t60;
t38 = 0.1e1 / t40 ^ 2;
t58 = t55 * t49;
t59 = t54 * t51;
t39 = t50 * t58 - t59;
t66 = t38 * t39;
t65 = t41 * t50;
t45 = t48 ^ 2;
t64 = t45 / t50 ^ 2;
t63 = t48 * t55;
t44 = 0.1e1 / (t54 ^ 2 * t64 + 0.1e1);
t62 = t54 * t44;
t56 = t39 ^ 2 * t38 + 0.1e1;
t46 = 0.1e1 / t50;
t37 = 0.1e1 / t40;
t36 = (0.1e1 + t64) * t62;
t35 = 0.1e1 / t56;
t32 = 0.1e1 / t34;
t31 = 0.1e1 / (t45 * t67 + 0.1e1);
t1 = [t46 * t44 * t63, 0, t36, 0, 0, 0; (-t32 * t61 - (-t42 * t45 * t46 * t62 + (t44 - 0.1e1) * t48 * t41) * t48 * t67) * t31, 0 (t50 * t32 - (-t54 * t65 + t42 * t48 + (-t42 * t61 + t65) * t36) * t48 * t33) * t55 * t31, 0, 0, 0; ((-t50 * t60 - t57) * t37 - (-t50 * t59 + t58) * t66) * t35, 0 (-t37 * t49 + t51 * t66) * t35 * t63, t56 * t35, 0, 0;];
Ja_rot  = t1;
