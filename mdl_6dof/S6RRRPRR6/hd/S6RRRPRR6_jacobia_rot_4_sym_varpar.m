% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR6_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:47
% EndTime: 2019-02-26 22:18:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (157->21), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->34)
t48 = cos(qJ(2));
t46 = sin(qJ(2));
t47 = sin(qJ(1));
t54 = t47 * t46;
t38 = atan2(-t54, -t48);
t36 = sin(t38);
t37 = cos(t38);
t30 = -t36 * t54 - t37 * t48;
t29 = 0.1e1 / t30 ^ 2;
t49 = cos(qJ(1));
t59 = t29 * t49 ^ 2;
t42 = qJ(3) + pkin(11);
t40 = sin(t42);
t41 = cos(t42);
t51 = t49 * t41;
t35 = t47 * t40 + t48 * t51;
t33 = 0.1e1 / t35 ^ 2;
t52 = t49 * t40;
t34 = -t47 * t41 + t48 * t52;
t58 = t33 * t34;
t43 = t46 ^ 2;
t57 = t43 / t48 ^ 2;
t56 = t46 * t49;
t39 = 0.1e1 / (t47 ^ 2 * t57 + 0.1e1);
t55 = t47 * t39;
t53 = t47 * t48;
t50 = t34 ^ 2 * t33 + 0.1e1;
t44 = 0.1e1 / t48;
t32 = 0.1e1 / t35;
t31 = (0.1e1 + t57) * t55;
t28 = 0.1e1 / t30;
t27 = 0.1e1 / t50;
t26 = 0.1e1 / (t43 * t59 + 0.1e1);
t1 = [t44 * t39 * t56, t31, 0, 0, 0, 0; (-t28 * t54 - (-t37 * t43 * t44 * t55 + (t39 - 0.1e1) * t46 * t36) * t46 * t59) * t26 (t48 * t28 - (-t36 * t53 + t37 * t46 + (t36 * t48 - t37 * t54) * t31) * t46 * t29) * t49 * t26, 0, 0, 0, 0; ((-t40 * t53 - t51) * t32 - (-t41 * t53 + t52) * t58) * t27 (-t32 * t40 + t41 * t58) * t27 * t56, t50 * t27, 0, 0, 0;];
Ja_rot  = t1;
