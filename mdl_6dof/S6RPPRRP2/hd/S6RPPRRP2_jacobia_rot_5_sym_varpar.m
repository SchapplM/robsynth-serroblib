% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:41
% EndTime: 2019-02-26 20:30:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (347->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t46 = pkin(10) + qJ(4);
t44 = cos(t46);
t42 = sin(t46);
t47 = qJ(1) + pkin(9);
t43 = sin(t47);
t55 = t43 * t42;
t37 = atan2(-t55, -t44);
t35 = sin(t37);
t36 = cos(t37);
t28 = -t35 * t55 - t36 * t44;
t27 = 0.1e1 / t28 ^ 2;
t45 = cos(t47);
t61 = t27 * t45 ^ 2;
t49 = cos(qJ(5));
t51 = t45 * t49;
t48 = sin(qJ(5));
t54 = t43 * t48;
t34 = t44 * t51 + t54;
t32 = 0.1e1 / t34 ^ 2;
t52 = t45 * t48;
t53 = t43 * t49;
t33 = t44 * t52 - t53;
t60 = t32 * t33;
t59 = t35 * t44;
t39 = t42 ^ 2;
t58 = t39 / t44 ^ 2;
t57 = t42 * t45;
t38 = 0.1e1 / (t43 ^ 2 * t58 + 0.1e1);
t56 = t43 * t38;
t50 = t33 ^ 2 * t32 + 0.1e1;
t40 = 0.1e1 / t44;
t31 = 0.1e1 / t34;
t30 = 0.1e1 / t50;
t29 = (0.1e1 + t58) * t56;
t26 = 0.1e1 / t28;
t25 = 0.1e1 / (t39 * t61 + 0.1e1);
t1 = [t40 * t38 * t57, 0, 0, t29, 0, 0; (-t26 * t55 - (-t36 * t39 * t40 * t56 + (t38 - 0.1e1) * t42 * t35) * t42 * t61) * t25, 0, 0 (t44 * t26 - (-t43 * t59 + t36 * t42 + (-t36 * t55 + t59) * t29) * t42 * t27) * t45 * t25, 0, 0; ((-t44 * t54 - t51) * t31 - (-t44 * t53 + t52) * t60) * t30, 0, 0 (-t31 * t48 + t49 * t60) * t30 * t57, t50 * t30, 0;];
Ja_rot  = t1;
