% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:01
% EndTime: 2019-02-26 21:02:01
% DurationCPUTime: 0.11s
% Computational Cost: add. (378->27), mult. (509->71), div. (100->11), fcn. (770->9), ass. (0->39)
t54 = sin(qJ(3));
t67 = t54 ^ 2;
t47 = qJ(1) + pkin(10);
t45 = sin(t47);
t46 = cos(t47);
t55 = cos(qJ(4));
t53 = sin(qJ(4));
t56 = cos(qJ(3));
t60 = t53 * t56;
t37 = t45 * t60 + t46 * t55;
t59 = t54 * t53;
t34 = atan2(-t37, t59);
t31 = sin(t34);
t32 = cos(t34);
t29 = -t31 * t37 + t32 * t59;
t28 = 0.1e1 / t29 ^ 2;
t40 = -t45 * t55 + t46 * t60;
t66 = t28 * t40;
t64 = t32 * t37;
t63 = t40 ^ 2 * t28;
t62 = t46 * t54;
t48 = 0.1e1 / t53;
t51 = 0.1e1 / t54;
t61 = t48 * t51;
t58 = t55 * t56;
t41 = t45 * t53 + t46 * t58;
t36 = 0.1e1 / t41 ^ 2;
t57 = t46 ^ 2 * t67 * t36;
t52 = 0.1e1 / t67;
t49 = 0.1e1 / t53 ^ 2;
t39 = t45 * t58 - t46 * t53;
t35 = 0.1e1 / t41;
t33 = 0.1e1 / (t37 ^ 2 * t52 * t49 + 0.1e1);
t30 = 0.1e1 / (0.1e1 + t57);
t27 = 0.1e1 / t29;
t26 = (t37 * t48 * t52 * t56 + t45) * t33;
t25 = 0.1e1 / (0.1e1 + t63);
t24 = (t37 * t49 * t55 - t39 * t48) * t51 * t33;
t1 = [-t40 * t33 * t61, 0, t26, t24, 0, 0; (-t37 * t27 - (-t31 + (t61 * t64 + t31) * t33) * t63) * t25, 0 (t26 * t64 * t66 + (-t27 * t62 - (t32 * t56 + (-t26 + t45) * t54 * t31) * t66) * t53) * t25 (t41 * t27 - (t32 * t54 * t55 - t31 * t39 + (-t31 * t59 - t64) * t24) * t66) * t25, 0, 0; (-t36 * t39 * t46 + t35 * t45) * t54 * t30, 0 (-t35 * t46 * t56 - t55 * t57) * t30, -t40 * t36 * t30 * t62, 0, 0;];
Ja_rot  = t1;
