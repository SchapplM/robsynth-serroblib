% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP7_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:19
% EndTime: 2019-02-26 21:11:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (238->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t43 = pkin(10) + qJ(3);
t42 = cos(t43);
t41 = sin(t43);
t45 = sin(qJ(1));
t53 = t45 * t41;
t36 = atan2(-t53, -t42);
t32 = sin(t36);
t33 = cos(t36);
t27 = -t32 * t53 - t33 * t42;
t26 = 0.1e1 / t27 ^ 2;
t47 = cos(qJ(1));
t59 = t26 * t47 ^ 2;
t46 = cos(qJ(4));
t49 = t47 * t46;
t44 = sin(qJ(4));
t52 = t45 * t44;
t35 = t42 * t49 + t52;
t31 = 0.1e1 / t35 ^ 2;
t50 = t47 * t44;
t51 = t45 * t46;
t34 = t42 * t50 - t51;
t58 = t31 * t34;
t57 = t32 * t42;
t38 = t41 ^ 2;
t56 = t38 / t42 ^ 2;
t55 = t41 * t47;
t37 = 0.1e1 / (t45 ^ 2 * t56 + 0.1e1);
t54 = t45 * t37;
t48 = t34 ^ 2 * t31 + 0.1e1;
t39 = 0.1e1 / t42;
t30 = 0.1e1 / t35;
t29 = 0.1e1 / t48;
t28 = (0.1e1 + t56) * t54;
t25 = 0.1e1 / t27;
t24 = 0.1e1 / (t38 * t59 + 0.1e1);
t1 = [t39 * t37 * t55, 0, t28, 0, 0, 0; (-t25 * t53 - (-t33 * t38 * t39 * t54 + (t37 - 0.1e1) * t41 * t32) * t41 * t59) * t24, 0 (t42 * t25 - (-t45 * t57 + t33 * t41 + (-t33 * t53 + t57) * t28) * t41 * t26) * t47 * t24, 0, 0, 0; ((-t42 * t52 - t49) * t30 - (-t42 * t51 + t50) * t58) * t29, 0 (-t30 * t44 + t46 * t58) * t29 * t55, t48 * t29, 0, 0;];
Ja_rot  = t1;
