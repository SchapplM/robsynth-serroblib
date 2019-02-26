% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:29
% EndTime: 2019-02-26 20:25:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (395->23), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->36)
t52 = pkin(10) + qJ(4);
t49 = cos(t52);
t46 = sin(t52);
t53 = qJ(1) + pkin(9);
t47 = sin(t53);
t58 = t47 * t46;
t40 = atan2(-t58, -t49);
t38 = sin(t40);
t39 = cos(t40);
t31 = -t38 * t58 - t39 * t49;
t30 = 0.1e1 / t31 ^ 2;
t50 = cos(t53);
t63 = t30 * t50 ^ 2;
t51 = pkin(11) + qJ(6);
t45 = sin(t51);
t48 = cos(t51);
t55 = t50 * t48;
t37 = t47 * t45 + t49 * t55;
t35 = 0.1e1 / t37 ^ 2;
t56 = t50 * t45;
t36 = -t47 * t48 + t49 * t56;
t62 = t35 * t36;
t42 = t46 ^ 2;
t61 = t42 / t49 ^ 2;
t60 = t46 * t50;
t41 = 0.1e1 / (t47 ^ 2 * t61 + 0.1e1);
t59 = t47 * t41;
t57 = t47 * t49;
t54 = t36 ^ 2 * t35 + 0.1e1;
t43 = 0.1e1 / t49;
t34 = 0.1e1 / t37;
t33 = (0.1e1 + t61) * t59;
t32 = 0.1e1 / t54;
t29 = 0.1e1 / t31;
t28 = 0.1e1 / (t42 * t63 + 0.1e1);
t1 = [t43 * t41 * t60, 0, 0, t33, 0, 0; (-t29 * t58 - (-t39 * t42 * t43 * t59 + (t41 - 0.1e1) * t46 * t38) * t46 * t63) * t28, 0, 0 (t49 * t29 - (-t38 * t57 + t39 * t46 + (t38 * t49 - t39 * t58) * t33) * t46 * t30) * t50 * t28, 0, 0; ((-t45 * t57 - t55) * t34 - (-t48 * t57 + t56) * t62) * t32, 0, 0 (-t34 * t45 + t48 * t62) * t32 * t60, 0, t54 * t32;];
Ja_rot  = t1;
