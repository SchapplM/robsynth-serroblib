% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:47
% EndTime: 2019-02-26 21:26:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (111->24), mult. (347->67), div. (52->9), fcn. (503->11), ass. (0->39)
t49 = cos(pkin(9));
t52 = sin(qJ(1));
t54 = cos(qJ(2));
t48 = sin(pkin(9));
t55 = cos(qJ(1));
t58 = t55 * t48;
t39 = -t52 * t49 + t54 * t58;
t57 = t55 * t49;
t40 = t52 * t48 + t54 * t57;
t50 = sin(qJ(5));
t53 = cos(qJ(5));
t32 = t39 * t50 + t40 * t53;
t30 = 0.1e1 / t32 ^ 2;
t31 = -t39 * t53 + t40 * t50;
t65 = t30 * t31;
t51 = sin(qJ(2));
t60 = t52 * t51;
t44 = atan2(t60, t54);
t41 = sin(t44);
t42 = cos(t44);
t35 = t41 * t60 + t42 * t54;
t34 = 0.1e1 / t35 ^ 2;
t64 = t34 * t55 ^ 2;
t45 = t51 ^ 2;
t63 = t45 / t54 ^ 2;
t62 = t51 * t55;
t43 = 0.1e1 / (t52 ^ 2 * t63 + 0.1e1);
t61 = t52 * t43;
t59 = t52 * t54;
t56 = t31 ^ 2 * t30 + 0.1e1;
t46 = 0.1e1 / t54;
t38 = -t49 * t59 + t58;
t37 = -t48 * t59 - t57;
t36 = (0.1e1 + t63) * t61;
t33 = 0.1e1 / t35;
t29 = 0.1e1 / t32;
t28 = 0.1e1 / (t45 * t64 + 0.1e1);
t27 = 0.1e1 / t56;
t1 = [t46 * t43 * t62, t36, 0, 0, 0, 0; (t33 * t60 + (t42 * t45 * t46 * t61 + (-t43 + 0.1e1) * t51 * t41) * t51 * t64) * t28 (-t54 * t33 + (t41 * t59 - t42 * t51 + (-t41 * t54 + t42 * t60) * t36) * t51 * t34) * t55 * t28, 0, 0, 0, 0; ((-t37 * t53 + t38 * t50) * t29 - (t37 * t50 + t38 * t53) * t65) * t27 ((t48 * t53 - t49 * t50) * t29 - (-t48 * t50 - t49 * t53) * t65) * t27 * t62, 0, 0, t56 * t27, 0;];
Ja_rot  = t1;
