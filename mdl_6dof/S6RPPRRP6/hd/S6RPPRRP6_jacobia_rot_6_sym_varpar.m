% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_jacobia_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:33:08
% EndTime: 2019-02-26 20:33:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (160->25), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->40)
t54 = cos(qJ(4));
t68 = t54 ^ 2;
t51 = sin(qJ(4));
t53 = cos(qJ(5));
t55 = cos(qJ(1));
t57 = t55 * t53;
t50 = sin(qJ(5));
t52 = sin(qJ(1));
t62 = t52 * t50;
t38 = t51 * t62 - t57;
t60 = t54 * t50;
t35 = atan2(-t38, t60);
t31 = sin(t35);
t32 = cos(t35);
t30 = -t31 * t38 + t32 * t60;
t29 = 0.1e1 / t30 ^ 2;
t58 = t55 * t50;
t61 = t52 * t53;
t41 = t51 * t58 + t61;
t67 = t29 * t41;
t65 = t32 * t38;
t64 = t41 ^ 2 * t29;
t44 = 0.1e1 / t50;
t47 = 0.1e1 / t54;
t63 = t44 * t47;
t59 = t54 * t55;
t42 = t51 * t57 - t62;
t37 = 0.1e1 / t42 ^ 2;
t56 = t55 ^ 2 * t68 * t37;
t48 = 0.1e1 / t68;
t45 = 0.1e1 / t50 ^ 2;
t40 = t51 * t61 + t58;
t36 = 0.1e1 / t42;
t34 = 0.1e1 / (t38 ^ 2 * t48 * t45 + 0.1e1);
t33 = 0.1e1 / (0.1e1 + t56);
t28 = 0.1e1 / t30;
t27 = (-t38 * t44 * t48 * t51 - t52) * t34;
t26 = 0.1e1 / (0.1e1 + t64);
t25 = (t38 * t45 * t53 - t40 * t44) * t47 * t34;
t1 = [-t41 * t34 * t63, 0, 0, t27, t25, 0; (-t38 * t28 - (-t31 + (t63 * t65 + t31) * t34) * t64) * t26, 0, 0 (t27 * t65 * t67 + (t28 * t59 - (-t32 * t51 + (-t27 - t52) * t31 * t54) * t67) * t50) * t26 (t42 * t28 - (t32 * t54 * t53 - t31 * t40 + (-t31 * t60 - t65) * t25) * t67) * t26, 0; (t37 * t40 * t55 - t36 * t52) * t54 * t33, 0, 0 (-t36 * t51 * t55 - t53 * t56) * t33, t41 * t37 * t33 * t59, 0;];
Ja_rot  = t1;
