% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:16
% EndTime: 2019-02-26 20:26:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (325->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t47 = pkin(10) + qJ(4);
t43 = sin(t47);
t48 = qJ(1) + pkin(9);
t44 = sin(t48);
t45 = cos(t47);
t57 = t44 * t45;
t38 = atan2(-t57, t43);
t36 = sin(t38);
t37 = cos(t38);
t30 = -t36 * t57 + t37 * t43;
t28 = 0.1e1 / t30 ^ 2;
t46 = cos(t48);
t62 = t28 * t46 ^ 2;
t49 = sin(qJ(6));
t53 = t46 * t49;
t50 = cos(qJ(6));
t55 = t44 * t50;
t35 = t43 * t53 + t55;
t33 = 0.1e1 / t35 ^ 2;
t52 = t46 * t50;
t56 = t44 * t49;
t34 = -t43 * t52 + t56;
t61 = t33 * t34;
t60 = t36 * t43;
t42 = t45 ^ 2;
t59 = 0.1e1 / t43 ^ 2 * t42;
t39 = 0.1e1 / (t44 ^ 2 * t59 + 0.1e1);
t58 = t44 * t39;
t54 = t45 * t46;
t51 = t34 ^ 2 * t33 + 0.1e1;
t40 = 0.1e1 / t43;
t32 = 0.1e1 / t35;
t31 = 0.1e1 / t51;
t29 = (0.1e1 + t59) * t58;
t27 = 0.1e1 / t30;
t26 = 0.1e1 / (t42 * t62 + 0.1e1);
t1 = [-t40 * t39 * t54, 0, 0, t29, 0, 0; (-t27 * t57 - (t37 * t40 * t42 * t58 + (t39 - 0.1e1) * t45 * t36) * t45 * t62) * t26, 0, 0 (-t43 * t27 - (t44 * t60 + t37 * t45 + (-t37 * t57 - t60) * t29) * t45 * t28) * t46 * t26, 0, 0; ((t43 * t55 + t53) * t32 - (-t43 * t56 + t52) * t61) * t31, 0, 0 (-t32 * t50 - t49 * t61) * t31 * t54, 0, t51 * t31;];
Ja_rot  = t1;
