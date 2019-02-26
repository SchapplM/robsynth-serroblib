% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:30
% EndTime: 2019-02-26 21:12:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (158->20), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->38)
t54 = sin(qJ(1));
t69 = t54 ^ 2;
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t56 = cos(qJ(1));
t58 = t56 * t55;
t45 = atan2(-t58, t53);
t43 = sin(t45);
t44 = cos(t45);
t37 = -t43 * t58 + t44 * t53;
t36 = 0.1e1 / t37 ^ 2;
t68 = t36 * t55;
t52 = qJ(4) + qJ(5);
t47 = sin(t52);
t60 = t56 * t47;
t48 = cos(t52);
t63 = t54 * t48;
t42 = t53 * t63 + t60;
t40 = 0.1e1 / t42 ^ 2;
t59 = t56 * t48;
t64 = t54 * t47;
t41 = t53 * t64 - t59;
t67 = t40 * t41;
t66 = t43 * t53;
t51 = t55 ^ 2;
t65 = 0.1e1 / t53 ^ 2 * t51;
t62 = t54 * t55;
t46 = 0.1e1 / (t56 ^ 2 * t65 + 0.1e1);
t61 = t56 * t46;
t57 = t41 ^ 2 * t40 + 0.1e1;
t49 = 0.1e1 / t53;
t39 = 0.1e1 / t42;
t38 = (0.1e1 + t65) * t61;
t35 = 0.1e1 / t37;
t34 = 0.1e1 / t57;
t33 = 0.1e1 / (t69 * t51 * t36 + 0.1e1);
t32 = t57 * t34;
t1 = [t49 * t46 * t62, 0, t38, 0, 0, 0; (-t35 * t58 + (-t44 * t49 * t51 * t61 + (-t46 + 0.1e1) * t55 * t43) * t69 * t68) * t33, 0 (t53 * t35 + (t56 * t66 + t44 * t55 + (-t44 * t58 - t66) * t38) * t68) * t54 * t33, 0, 0, 0; ((t53 * t60 + t63) * t39 - (t53 * t59 - t64) * t67) * t34, 0 (t39 * t47 - t48 * t67) * t34 * t62, t32, t32, 0;];
Ja_rot  = t1;
