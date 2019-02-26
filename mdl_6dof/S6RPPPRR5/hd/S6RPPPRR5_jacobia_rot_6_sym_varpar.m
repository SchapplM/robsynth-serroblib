% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPPRR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:53
% EndTime: 2019-02-26 20:24:53
% DurationCPUTime: 0.09s
% Computational Cost: add. (218->21), mult. (442->61), div. (52->9), fcn. (659->11), ass. (0->38)
t45 = cos(pkin(9));
t50 = cos(qJ(1));
t52 = sin(pkin(9));
t62 = sin(qJ(1));
t39 = t62 * t45 + t50 * t52;
t63 = t39 ^ 2;
t49 = cos(qJ(5));
t38 = t50 * t45 - t62 * t52;
t47 = sin(qJ(5));
t57 = t38 * t47;
t36 = atan2(t57, -t49);
t34 = sin(t36);
t35 = cos(t36);
t29 = t34 * t57 - t35 * t49;
t28 = 0.1e1 / t29 ^ 2;
t61 = t28 * t47;
t46 = sin(qJ(6));
t48 = cos(qJ(6));
t53 = t48 * t49;
t33 = -t38 * t46 + t39 * t53;
t31 = 0.1e1 / t33 ^ 2;
t54 = t46 * t49;
t32 = t38 * t48 + t39 * t54;
t60 = t31 * t32;
t59 = t34 * t49;
t42 = t47 ^ 2;
t55 = t42 / t49 ^ 2;
t37 = 0.1e1 / (t38 ^ 2 * t55 + 0.1e1);
t58 = t38 * t37;
t56 = t39 * t47;
t51 = t32 ^ 2 * t31 + 0.1e1;
t43 = 0.1e1 / t49;
t30 = 0.1e1 / t33;
t27 = 0.1e1 / t29;
t26 = (-0.1e1 - t55) * t58;
t25 = 0.1e1 / t51;
t24 = 0.1e1 / (t63 * t42 * t28 + 0.1e1);
t1 = [t43 * t37 * t56, 0, 0, 0, t26, 0; (t27 * t57 - (t35 * t42 * t43 * t58 + (t37 - 0.1e1) * t47 * t34) * t63 * t61) * t24, 0, 0, 0 (t49 * t27 - (t38 * t59 + t35 * t47 + (t35 * t57 + t59) * t26) * t61) * t39 * t24, 0; ((t38 * t54 - t39 * t48) * t30 - (t38 * t53 + t39 * t46) * t60) * t25, 0, 0, 0 (-t30 * t46 + t48 * t60) * t25 * t56, t51 * t25;];
Ja_rot  = t1;
