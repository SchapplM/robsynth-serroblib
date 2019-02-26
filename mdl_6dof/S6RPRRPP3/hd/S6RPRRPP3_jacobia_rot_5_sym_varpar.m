% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RPRRPP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:57:35
% EndTime: 2019-02-26 20:57:35
% DurationCPUTime: 0.09s
% Computational Cost: add. (369->27), mult. (486->71), div. (108->11), fcn. (755->9), ass. (0->38)
t46 = qJ(1) + pkin(9);
t44 = sin(t46);
t45 = cos(t46);
t53 = cos(qJ(4));
t51 = sin(qJ(4));
t54 = cos(qJ(3));
t57 = t51 * t54;
t35 = t44 * t57 + t45 * t53;
t52 = sin(qJ(3));
t56 = t52 * t51;
t34 = atan2(-t35, t56);
t31 = sin(t34);
t32 = cos(t34);
t29 = -t31 * t35 + t32 * t56;
t28 = 0.1e1 / t29 ^ 2;
t38 = -t44 * t53 + t45 * t57;
t64 = t28 * t38;
t55 = t53 * t54;
t39 = t44 * t51 + t45 * t55;
t43 = 0.1e1 / t45 ^ 2;
t50 = 0.1e1 / t52 ^ 2;
t30 = 0.1e1 / (t39 ^ 2 * t43 * t50 + 0.1e1);
t49 = 0.1e1 / t52;
t63 = t30 * t49;
t61 = t32 * t35;
t60 = t38 ^ 2 * t28;
t47 = 0.1e1 / t51;
t59 = t47 * t49;
t58 = t50 * t54;
t48 = 0.1e1 / t51 ^ 2;
t42 = 0.1e1 / t45;
t37 = t44 * t55 - t45 * t51;
t33 = 0.1e1 / (t35 ^ 2 * t50 * t48 + 0.1e1);
t27 = 0.1e1 / t29;
t26 = (t35 * t47 * t58 + t44) * t33;
t25 = 0.1e1 / (0.1e1 + t60);
t24 = (t35 * t48 * t53 - t37 * t47) * t49 * t33;
t1 = [-t38 * t33 * t59, 0, t26, t24, 0, 0; (-t35 * t27 - (-t31 + (t59 * t61 + t31) * t33) * t60) * t25, 0 (t26 * t61 * t64 + (-t45 * t52 * t27 - (t32 * t54 + (-t26 + t44) * t52 * t31) * t64) * t51) * t25 (t39 * t27 - (t32 * t52 * t53 - t31 * t37 + (-t31 * t56 - t61) * t24) * t64) * t25, 0, 0; (t39 * t43 * t44 - t37 * t42) * t63, 0 (-t39 * t42 * t58 - t53) * t30, -t38 * t42 * t63, 0, 0;];
Ja_rot  = t1;
