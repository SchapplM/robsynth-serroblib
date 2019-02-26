% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:57:25
% EndTime: 2019-02-26 20:57:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (369->27), mult. (486->71), div. (108->11), fcn. (755->9), ass. (0->38)
t45 = qJ(1) + pkin(9);
t43 = sin(t45);
t44 = cos(t45);
t50 = sin(qJ(4));
t52 = cos(qJ(4));
t53 = cos(qJ(3));
t54 = t52 * t53;
t36 = t43 * t54 - t44 * t50;
t51 = sin(qJ(3));
t55 = t51 * t52;
t34 = atan2(-t36, t55);
t31 = sin(t34);
t32 = cos(t34);
t29 = -t31 * t36 + t32 * t55;
t28 = 0.1e1 / t29 ^ 2;
t39 = t43 * t50 + t44 * t54;
t63 = t28 * t39;
t56 = t50 * t53;
t38 = t43 * t52 - t44 * t56;
t42 = 0.1e1 / t44 ^ 2;
t47 = 0.1e1 / t51 ^ 2;
t30 = 0.1e1 / (t38 ^ 2 * t42 * t47 + 0.1e1);
t46 = 0.1e1 / t51;
t62 = t30 * t46;
t60 = t32 * t36;
t59 = t39 ^ 2 * t28;
t48 = 0.1e1 / t52;
t58 = t46 * t48;
t57 = t47 * t53;
t49 = 0.1e1 / t52 ^ 2;
t41 = 0.1e1 / t44;
t35 = t43 * t56 + t44 * t52;
t33 = 0.1e1 / (t36 ^ 2 * t47 * t49 + 0.1e1);
t27 = 0.1e1 / t29;
t26 = (t36 * t48 * t57 + t43) * t33;
t25 = 0.1e1 / (0.1e1 + t59);
t24 = (-t36 * t49 * t50 + t35 * t48) * t46 * t33;
t1 = [-t39 * t33 * t58, 0, t26, t24, 0, 0; (-t36 * t27 - (-t31 + (t58 * t60 + t31) * t33) * t59) * t25, 0 (t26 * t60 * t63 + (-t44 * t51 * t27 - (t32 * t53 + (-t26 + t43) * t31 * t51) * t63) * t52) * t25 (t38 * t27 - (-t32 * t51 * t50 + t31 * t35 + (-t31 * t55 - t60) * t24) * t63) * t25, 0, 0; (t38 * t42 * t43 + t35 * t41) * t62, 0 (-t38 * t41 * t57 + t50) * t30, -t39 * t41 * t62, 0, 0;];
Ja_rot  = t1;
