% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:27
% EndTime: 2019-02-26 20:46:27
% DurationCPUTime: 0.10s
% Computational Cost: add. (216->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t44 = pkin(9) + qJ(3);
t42 = sin(t44);
t43 = cos(t44);
t46 = sin(qJ(1));
t54 = t46 * t43;
t37 = atan2(-t54, t42);
t35 = sin(t37);
t36 = cos(t37);
t28 = -t35 * t54 + t36 * t42;
t27 = 0.1e1 / t28 ^ 2;
t48 = cos(qJ(1));
t60 = t27 * t48 ^ 2;
t45 = sin(qJ(5));
t51 = t48 * t45;
t47 = cos(qJ(5));
t52 = t46 * t47;
t34 = t42 * t51 + t52;
t32 = 0.1e1 / t34 ^ 2;
t50 = t48 * t47;
t53 = t46 * t45;
t33 = -t42 * t50 + t53;
t59 = t32 * t33;
t58 = t35 * t42;
t41 = t43 ^ 2;
t57 = 0.1e1 / t42 ^ 2 * t41;
t56 = t43 * t48;
t38 = 0.1e1 / (t46 ^ 2 * t57 + 0.1e1);
t55 = t46 * t38;
t49 = t33 ^ 2 * t32 + 0.1e1;
t39 = 0.1e1 / t42;
t31 = 0.1e1 / t34;
t30 = 0.1e1 / t49;
t29 = (0.1e1 + t57) * t55;
t26 = 0.1e1 / t28;
t25 = 0.1e1 / (t41 * t60 + 0.1e1);
t1 = [-t39 * t38 * t56, 0, t29, 0, 0, 0; (-t26 * t54 - (t36 * t39 * t41 * t55 + (t38 - 0.1e1) * t43 * t35) * t43 * t60) * t25, 0 (-t42 * t26 - (t46 * t58 + t36 * t43 + (-t36 * t54 - t58) * t29) * t43 * t27) * t48 * t25, 0, 0, 0; ((t42 * t52 + t51) * t31 - (-t42 * t53 + t50) * t59) * t30, 0 (-t31 * t47 - t45 * t59) * t30 * t56, 0, t49 * t30, 0;];
Ja_rot  = t1;
