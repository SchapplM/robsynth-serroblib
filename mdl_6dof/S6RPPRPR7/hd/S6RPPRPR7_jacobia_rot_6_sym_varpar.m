% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:09
% EndTime: 2019-02-26 20:29:09
% DurationCPUTime: 0.10s
% Computational Cost: add. (263->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t50 = pkin(9) + qJ(4);
t46 = sin(t50);
t48 = cos(t50);
t52 = cos(qJ(1));
t54 = t52 * t48;
t40 = atan2(-t54, t46);
t38 = sin(t40);
t39 = cos(t40);
t31 = -t38 * t54 + t39 * t46;
t30 = 0.1e1 / t31 ^ 2;
t51 = sin(qJ(1));
t64 = t30 * t51 ^ 2;
t49 = pkin(10) + qJ(6);
t45 = sin(t49);
t56 = t52 * t45;
t47 = cos(t49);
t58 = t51 * t47;
t37 = t46 * t58 + t56;
t35 = 0.1e1 / t37 ^ 2;
t55 = t52 * t47;
t59 = t51 * t45;
t36 = t46 * t59 - t55;
t63 = t35 * t36;
t62 = t38 * t46;
t44 = t48 ^ 2;
t61 = 0.1e1 / t46 ^ 2 * t44;
t60 = t48 * t51;
t41 = 0.1e1 / (t52 ^ 2 * t61 + 0.1e1);
t57 = t52 * t41;
t53 = t36 ^ 2 * t35 + 0.1e1;
t42 = 0.1e1 / t46;
t34 = 0.1e1 / t37;
t33 = (0.1e1 + t61) * t57;
t32 = 0.1e1 / t53;
t29 = 0.1e1 / t31;
t28 = 0.1e1 / (t44 * t64 + 0.1e1);
t1 = [t42 * t41 * t60, 0, 0, t33, 0, 0; (-t29 * t54 + (-t39 * t42 * t44 * t57 + (-t41 + 0.1e1) * t48 * t38) * t48 * t64) * t28, 0, 0 (t46 * t29 + (t52 * t62 + t39 * t48 + (-t39 * t54 - t62) * t33) * t48 * t30) * t51 * t28, 0, 0; ((t46 * t56 + t58) * t34 - (t46 * t55 - t59) * t63) * t32, 0, 0 (t34 * t45 - t47 * t63) * t32 * t60, 0, t53 * t32;];
Ja_rot  = t1;
