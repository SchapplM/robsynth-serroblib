% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPPRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:17
% EndTime: 2019-02-26 20:24:17
% DurationCPUTime: 0.09s
% Computational Cost: add. (196->23), mult. (442->61), div. (52->9), fcn. (659->11), ass. (0->38)
t53 = sin(pkin(9));
t54 = cos(pkin(9));
t64 = sin(qJ(1));
t65 = cos(qJ(1));
t41 = t65 * t53 - t64 * t54;
t66 = t41 ^ 2;
t49 = sin(qJ(5));
t40 = -t64 * t53 - t65 * t54;
t51 = cos(qJ(5));
t59 = t40 * t51;
t38 = atan2(-t59, -t49);
t36 = sin(t38);
t37 = cos(t38);
t31 = -t36 * t59 - t37 * t49;
t30 = 0.1e1 / t31 ^ 2;
t63 = t30 * t51;
t48 = sin(qJ(6));
t50 = cos(qJ(6));
t55 = t49 * t50;
t35 = -t40 * t48 + t41 * t55;
t33 = 0.1e1 / t35 ^ 2;
t56 = t48 * t49;
t34 = t40 * t50 + t41 * t56;
t62 = t33 * t34;
t61 = t36 * t49;
t47 = t51 ^ 2;
t57 = 0.1e1 / t49 ^ 2 * t47;
t39 = 0.1e1 / (t40 ^ 2 * t57 + 0.1e1);
t60 = t40 * t39;
t58 = t41 * t51;
t52 = t34 ^ 2 * t33 + 0.1e1;
t45 = 0.1e1 / t49;
t32 = 0.1e1 / t35;
t29 = 0.1e1 / t31;
t28 = (-0.1e1 - t57) * t60;
t27 = 0.1e1 / t52;
t26 = 0.1e1 / (t66 * t47 * t30 + 0.1e1);
t1 = [-t45 * t39 * t58, 0, 0, 0, t28, 0; (-t29 * t59 + (t37 * t45 * t47 * t60 + (-t39 + 0.1e1) * t51 * t36) * t66 * t63) * t26, 0, 0, 0 (t49 * t29 + (t40 * t61 - t37 * t51 + (-t37 * t59 + t61) * t28) * t63) * t41 * t26, 0; ((t40 * t56 - t41 * t50) * t32 - (t40 * t55 + t41 * t48) * t62) * t27, 0, 0, 0 (t32 * t48 - t50 * t62) * t27 * t58, t52 * t27;];
Ja_rot  = t1;
