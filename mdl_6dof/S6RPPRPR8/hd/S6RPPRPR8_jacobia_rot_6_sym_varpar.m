% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:44
% EndTime: 2019-02-26 20:29:44
% DurationCPUTime: 0.07s
% Computational Cost: add. (193->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t43 = pkin(9) + qJ(4);
t42 = cos(t43);
t41 = sin(t43);
t47 = cos(qJ(1));
t51 = t47 * t41;
t37 = atan2(t51, t42);
t34 = sin(t37);
t35 = cos(t37);
t27 = t34 * t51 + t35 * t42;
t26 = 0.1e1 / t27 ^ 2;
t45 = sin(qJ(1));
t59 = t26 * t45 ^ 2;
t46 = cos(qJ(6));
t49 = t47 * t46;
t44 = sin(qJ(6));
t54 = t45 * t44;
t33 = -t42 * t54 + t49;
t31 = 0.1e1 / t33 ^ 2;
t50 = t47 * t44;
t53 = t45 * t46;
t32 = t42 * t53 + t50;
t58 = t31 * t32;
t57 = t34 * t42;
t38 = t41 ^ 2;
t56 = t38 / t42 ^ 2;
t55 = t41 * t45;
t36 = 0.1e1 / (t47 ^ 2 * t56 + 0.1e1);
t52 = t47 * t36;
t48 = t32 ^ 2 * t31 + 0.1e1;
t39 = 0.1e1 / t42;
t30 = 0.1e1 / t33;
t29 = 0.1e1 / t48;
t28 = (0.1e1 + t56) * t52;
t25 = 0.1e1 / t27;
t24 = 0.1e1 / (t38 * t59 + 0.1e1);
t1 = [-t39 * t36 * t55, 0, 0, t28, 0, 0; (t25 * t51 - (-t35 * t38 * t39 * t52 + (t36 - 0.1e1) * t41 * t34) * t41 * t59) * t24, 0, 0 (t42 * t25 - (t47 * t57 - t35 * t41 + (t35 * t51 - t57) * t28) * t41 * t26) * t45 * t24, 0, 0; ((t42 * t49 - t54) * t30 - (-t42 * t50 - t53) * t58) * t29, 0, 0 (-t30 * t46 - t44 * t58) * t29 * t55, 0, t48 * t29;];
Ja_rot  = t1;
