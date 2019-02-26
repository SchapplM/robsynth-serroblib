% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:17
% EndTime: 2019-02-26 20:45:18
% DurationCPUTime: 0.14s
% Computational Cost: add. (378->26), mult. (509->71), div. (100->11), fcn. (770->9), ass. (0->39)
t55 = cos(qJ(3));
t66 = t55 ^ 2;
t46 = qJ(1) + pkin(9);
t44 = sin(t46);
t45 = cos(t46);
t52 = sin(qJ(5));
t53 = sin(qJ(3));
t54 = cos(qJ(5));
t58 = t53 * t54;
t40 = t44 * t58 + t45 * t52;
t57 = t55 * t54;
t35 = atan2(t40, t57);
t32 = sin(t35);
t33 = cos(t35);
t30 = t32 * t40 + t33 * t57;
t29 = 0.1e1 / t30 ^ 2;
t38 = t44 * t52 - t45 * t58;
t65 = t29 * t38;
t63 = t33 * t40;
t62 = t38 ^ 2 * t29;
t61 = t45 * t55;
t47 = 0.1e1 / t54;
t50 = 0.1e1 / t55;
t60 = t47 * t50;
t59 = t52 * t53;
t39 = t44 * t54 + t45 * t59;
t37 = 0.1e1 / t39 ^ 2;
t56 = t45 ^ 2 * t66 * t37;
t51 = 0.1e1 / t66;
t48 = 0.1e1 / t54 ^ 2;
t41 = -t44 * t59 + t45 * t54;
t36 = 0.1e1 / t39;
t34 = 0.1e1 / (t40 ^ 2 * t51 * t48 + 0.1e1);
t31 = 0.1e1 / (0.1e1 + t56);
t28 = 0.1e1 / t30;
t27 = (t40 * t47 * t51 * t53 + t44) * t34;
t26 = 0.1e1 / (0.1e1 + t62);
t25 = (t40 * t48 * t52 + t41 * t47) * t50 * t34;
t1 = [-t38 * t34 * t60, 0, t27, 0, t25, 0; (t40 * t28 - (-t32 + (-t60 * t63 + t32) * t34) * t62) * t26, 0 (-t27 * t63 * t65 + (-t28 * t61 - (-t33 * t53 + (-t27 + t44) * t32 * t55) * t65) * t54) * t26, 0 (t39 * t28 - (-t33 * t55 * t52 + t32 * t41 + (-t32 * t57 + t63) * t25) * t65) * t26, 0; (t37 * t41 * t45 + t36 * t44) * t55 * t31, 0 (t36 * t45 * t53 + t52 * t56) * t31, 0, -t38 * t37 * t31 * t61, 0;];
Ja_rot  = t1;
