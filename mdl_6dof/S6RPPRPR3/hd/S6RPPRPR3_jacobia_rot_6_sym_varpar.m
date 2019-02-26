% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:47
% EndTime: 2019-02-26 20:26:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (324->21), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->38)
t50 = qJ(1) + pkin(9);
t46 = sin(t50);
t65 = t46 ^ 2;
t49 = qJ(4) + pkin(10);
t45 = sin(t49);
t47 = cos(t49);
t48 = cos(t50);
t56 = t48 * t47;
t40 = atan2(-t56, t45);
t38 = sin(t40);
t39 = cos(t40);
t31 = -t38 * t56 + t39 * t45;
t30 = 0.1e1 / t31 ^ 2;
t64 = t30 * t47;
t51 = sin(qJ(6));
t55 = t48 * t51;
t52 = cos(qJ(6));
t58 = t46 * t52;
t37 = t45 * t58 + t55;
t35 = 0.1e1 / t37 ^ 2;
t54 = t48 * t52;
t59 = t46 * t51;
t36 = t45 * t59 - t54;
t63 = t35 * t36;
t62 = t38 * t45;
t44 = t47 ^ 2;
t61 = 0.1e1 / t45 ^ 2 * t44;
t60 = t46 * t47;
t41 = 0.1e1 / (t48 ^ 2 * t61 + 0.1e1);
t57 = t48 * t41;
t53 = t36 ^ 2 * t35 + 0.1e1;
t42 = 0.1e1 / t45;
t34 = 0.1e1 / t37;
t33 = 0.1e1 / t53;
t32 = (0.1e1 + t61) * t57;
t29 = 0.1e1 / t31;
t28 = 0.1e1 / (t65 * t44 * t30 + 0.1e1);
t1 = [t42 * t41 * t60, 0, 0, t32, 0, 0; (-t29 * t56 + (-t39 * t42 * t44 * t57 + (-t41 + 0.1e1) * t47 * t38) * t65 * t64) * t28, 0, 0 (t45 * t29 + (t48 * t62 + t39 * t47 + (-t39 * t56 - t62) * t32) * t64) * t46 * t28, 0, 0; ((t45 * t55 + t58) * t34 - (t45 * t54 - t59) * t63) * t33, 0, 0 (t34 * t51 - t52 * t63) * t33 * t60, 0, t53 * t33;];
Ja_rot  = t1;
