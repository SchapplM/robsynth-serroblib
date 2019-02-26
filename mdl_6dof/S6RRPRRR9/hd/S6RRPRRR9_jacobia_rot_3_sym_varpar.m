% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR9_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobia_rot_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:45
% EndTime: 2019-02-26 21:58:45
% DurationCPUTime: 0.08s
% Computational Cost: add. (156->24), mult. (403->63), div. (67->11), fcn. (610->11), ass. (0->39)
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t49 = cos(qJ(1));
t47 = sin(qJ(1));
t52 = cos(pkin(6));
t51 = t47 * t52;
t38 = -t46 * t51 + t49 * t48;
t43 = sin(pkin(12));
t45 = cos(pkin(12));
t44 = sin(pkin(6));
t55 = t44 * t47;
t29 = t38 * t45 + t43 * t55;
t27 = 0.1e1 / t29 ^ 2;
t28 = t38 * t43 - t45 * t55;
t59 = t27 * t28;
t50 = t49 * t52;
t34 = t47 * t46 - t48 * t50;
t54 = t44 * t48;
t32 = atan2(-t34, -t54);
t31 = cos(t32);
t58 = t31 * t34;
t30 = sin(t32);
t24 = -t30 * t34 - t31 * t54;
t23 = 0.1e1 / t24 ^ 2;
t37 = t49 * t46 + t48 * t51;
t57 = t37 ^ 2 * t23;
t40 = 0.1e1 / t44;
t41 = 0.1e1 / t48;
t56 = t40 * t41;
t53 = t44 * t49;
t42 = 0.1e1 / t48 ^ 2;
t36 = t46 * t50 + t47 * t48;
t33 = 0.1e1 / (0.1e1 + t34 ^ 2 / t44 ^ 2 * t42);
t26 = 0.1e1 / t29;
t25 = 0.1e1 / (t28 ^ 2 * t27 + 0.1e1);
t22 = 0.1e1 / t24;
t21 = 0.1e1 / (0.1e1 + t57);
t20 = (t34 * t42 * t46 + t36 * t41) * t40 * t33;
t1 = [t37 * t33 * t56, t20, 0, 0, 0, 0; (-t34 * t22 - (-t30 + (-t56 * t58 + t30) * t33) * t57) * t21 (t38 * t22 - (t31 * t44 * t46 - t30 * t36 + (t30 * t54 - t58) * t20) * t37 * t23) * t21, 0, 0, 0, 0; ((-t36 * t43 - t45 * t53) * t26 - (-t36 * t45 + t43 * t53) * t59) * t25 (-t43 * t26 + t45 * t59) * t37 * t25, 0, 0, 0, 0;];
Ja_rot  = t1;
