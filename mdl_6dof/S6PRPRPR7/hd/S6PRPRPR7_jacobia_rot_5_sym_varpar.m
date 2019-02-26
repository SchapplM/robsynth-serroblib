% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:47
% EndTime: 2019-02-26 19:49:47
% DurationCPUTime: 0.09s
% Computational Cost: add. (251->25), mult. (724->66), div. (56->9), fcn. (1043->11), ass. (0->42)
t53 = sin(pkin(10));
t55 = cos(pkin(10));
t58 = sin(qJ(2));
t56 = cos(pkin(6));
t60 = cos(qJ(2));
t62 = t56 * t60;
t47 = t53 * t58 - t55 * t62;
t59 = cos(qJ(4));
t54 = sin(pkin(6));
t57 = sin(qJ(4));
t67 = t54 * t57;
t42 = t47 * t59 + t55 * t67;
t64 = t54 * t60;
t51 = t56 * t57 + t59 * t64;
t39 = atan2(t42, t51);
t35 = sin(t39);
t36 = cos(t39);
t34 = t35 * t42 + t36 * t51;
t33 = 0.1e1 / t34 ^ 2;
t49 = t53 * t62 + t55 * t58;
t40 = -t49 * t59 + t53 * t67;
t69 = t33 * t40;
t46 = 0.1e1 / t51 ^ 2;
t68 = t42 * t46;
t66 = t54 * t58;
t65 = t54 * t59;
t63 = t56 * t58;
t61 = -t35 * t51 + t36 * t42;
t52 = t56 * t59 - t57 * t64;
t50 = -t53 * t63 + t55 * t60;
t48 = t53 * t60 + t55 * t63;
t45 = 0.1e1 / t51;
t44 = 0.1e1 / t50 ^ 2;
t43 = -t47 * t57 + t55 * t65;
t41 = t49 * t57 + t53 * t65;
t38 = 0.1e1 / (t42 ^ 2 * t46 + 0.1e1);
t37 = 0.1e1 / (t41 ^ 2 * t44 + 0.1e1);
t32 = 0.1e1 / t34;
t31 = 0.1e1 / (t40 ^ 2 * t33 + 0.1e1);
t30 = (t45 * t48 + t66 * t68) * t59 * t38;
t29 = (t43 * t45 - t52 * t68) * t38;
t1 = [0, t30, 0, t29, 0, 0; 0 (-t50 * t59 * t32 - ((t35 * t48 - t36 * t66) * t59 + t61 * t30) * t69) * t31, 0 (t41 * t32 - (t61 * t29 + t35 * t43 + t36 * t52) * t69) * t31, 0, 0; 0 (t41 * t44 * t49 + t57) * t37, 0, -t40 / t50 * t37, 0, 0;];
Ja_rot  = t1;
