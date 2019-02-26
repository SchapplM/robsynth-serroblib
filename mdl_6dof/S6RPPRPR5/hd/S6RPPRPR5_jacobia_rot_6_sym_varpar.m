% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:00
% EndTime: 2019-02-26 20:28:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (111->19), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t44 = sin(qJ(4));
t45 = sin(qJ(1));
t46 = cos(qJ(4));
t52 = t45 * t46;
t37 = atan2(t52, t44);
t34 = sin(t37);
t35 = cos(t37);
t28 = t34 * t52 + t35 * t44;
t27 = 0.1e1 / t28 ^ 2;
t47 = cos(qJ(1));
t59 = t27 * t47 ^ 2;
t40 = pkin(9) + qJ(6);
t39 = cos(t40);
t49 = t47 * t39;
t38 = sin(t40);
t54 = t45 * t38;
t33 = t44 * t49 - t54;
t31 = 0.1e1 / t33 ^ 2;
t50 = t47 * t38;
t53 = t45 * t39;
t32 = t44 * t50 + t53;
t58 = t31 * t32;
t57 = t34 * t44;
t43 = t46 ^ 2;
t56 = 0.1e1 / t44 ^ 2 * t43;
t36 = 0.1e1 / (t45 ^ 2 * t56 + 0.1e1);
t55 = t45 * t36;
t51 = t46 * t47;
t48 = t32 ^ 2 * t31 + 0.1e1;
t41 = 0.1e1 / t44;
t30 = 0.1e1 / t33;
t29 = (-0.1e1 - t56) * t55;
t26 = 0.1e1 / t28;
t25 = 0.1e1 / (t43 * t59 + 0.1e1);
t24 = 0.1e1 / t48;
t1 = [t41 * t36 * t51, 0, 0, t29, 0, 0; (t26 * t52 + (t35 * t41 * t43 * t55 + (-t36 + 0.1e1) * t46 * t34) * t46 * t59) * t25, 0, 0 (t44 * t26 + (-t45 * t57 + t35 * t46 + (t35 * t52 - t57) * t29) * t46 * t27) * t47 * t25, 0, 0; ((-t44 * t54 + t49) * t30 - (-t44 * t53 - t50) * t58) * t24, 0, 0 (t30 * t38 - t39 * t58) * t24 * t51, 0, t48 * t24;];
Ja_rot  = t1;
