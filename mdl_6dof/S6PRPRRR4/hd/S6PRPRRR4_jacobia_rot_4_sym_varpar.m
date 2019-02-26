% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR4_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:55:33
% EndTime: 2019-02-26 19:55:33
% DurationCPUTime: 0.06s
% Computational Cost: add. (131->19), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->30)
t47 = sin(pkin(11));
t48 = sin(pkin(6));
t57 = t47 * t48;
t52 = cos(qJ(2));
t56 = t48 * t52;
t50 = cos(pkin(6));
t51 = sin(qJ(2));
t55 = t50 * t51;
t54 = t50 * t52;
t49 = cos(pkin(11));
t40 = -t47 * t55 + t49 * t52;
t45 = pkin(12) + qJ(4);
t42 = sin(t45);
t43 = cos(t45);
t31 = t40 * t43 + t42 * t57;
t29 = 0.1e1 / t31 ^ 2;
t30 = t40 * t42 - t43 * t57;
t53 = t30 ^ 2 * t29 + 0.1e1;
t46 = 0.1e1 / t52 ^ 2;
t39 = t47 * t54 + t49 * t51;
t38 = t47 * t52 + t49 * t55;
t36 = t47 * t51 - t49 * t54;
t34 = atan2(-t36, -t56);
t33 = cos(t34);
t32 = sin(t34);
t28 = 0.1e1 / t53;
t27 = -t32 * t36 - t33 * t56;
t26 = 0.1e1 / t27 ^ 2;
t24 = (t38 / t52 + t51 * t36 * t46) / t48 / (0.1e1 + t36 ^ 2 / t48 ^ 2 * t46);
t1 = [0, t24, 0, 0, 0, 0; 0 (t40 / t27 - (t33 * t48 * t51 - t32 * t38 + (t32 * t56 - t33 * t36) * t24) * t39 * t26) / (t39 ^ 2 * t26 + 0.1e1) 0, 0, 0, 0; 0 (-t42 / t31 + t43 * t30 * t29) * t39 * t28, 0, t53 * t28, 0, 0;];
Ja_rot  = t1;
