% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR7_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_rot_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:29
% EndTime: 2019-02-26 19:57:29
% DurationCPUTime: 0.08s
% Computational Cost: add. (148->22), mult. (435->55), div. (30->9), fcn. (617->13), ass. (0->34)
t43 = sin(pkin(13));
t47 = cos(pkin(13));
t51 = cos(qJ(2));
t49 = cos(pkin(6));
t50 = sin(qJ(2));
t54 = t49 * t50;
t40 = -t43 * t54 + t47 * t51;
t48 = cos(pkin(7));
t58 = t40 * t48;
t44 = sin(pkin(7));
t45 = sin(pkin(6));
t57 = t44 * t45;
t56 = t45 * t48;
t55 = t45 * t50;
t53 = t49 * t51;
t39 = -t43 * t53 - t47 * t50;
t52 = t39 * t48 + t43 * t57;
t46 = cos(pkin(14));
t42 = sin(pkin(14));
t38 = -t43 * t51 - t47 * t54;
t37 = t49 * t48 - t51 * t57;
t36 = 0.1e1 / t37 ^ 2;
t35 = -t39 * t44 + t43 * t56;
t34 = (-t43 * t50 + t47 * t53) * t44 + t47 * t56;
t33 = atan2(t34, t37);
t31 = cos(t33);
t30 = sin(t33);
t29 = t40 * t46 + t52 * t42;
t28 = t40 * t42 - t52 * t46;
t27 = 0.1e1 / t29 ^ 2;
t25 = t30 * t34 + t31 * t37;
t24 = 0.1e1 / t25 ^ 2;
t22 = (t38 / t37 - t34 * t36 * t55) * t44 / (t34 ^ 2 * t36 + 0.1e1);
t1 = [0, t22, 0, 0, 0, 0; 0 (t40 * t44 / t25 - ((t30 * t38 + t31 * t55) * t44 + (-t30 * t37 + t31 * t34) * t22) * t35 * t24) / (t35 ^ 2 * t24 + 0.1e1) 0, 0, 0, 0; 0 ((t39 * t42 + t46 * t58) / t29 - (t39 * t46 - t42 * t58) * t28 * t27) / (t28 ^ 2 * t27 + 0.1e1) 0, 0, 0, 0;];
Ja_rot  = t1;
