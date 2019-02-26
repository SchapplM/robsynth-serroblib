% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR8_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:34
% EndTime: 2019-02-26 21:32:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (102->20), mult. (329->52), div. (61->11), fcn. (494->9), ass. (0->34)
t45 = sin(qJ(2));
t59 = t45 ^ 2;
t43 = sin(pkin(10));
t44 = cos(pkin(10));
t48 = cos(qJ(1));
t50 = t48 * t44;
t46 = sin(qJ(1));
t47 = cos(qJ(2));
t52 = t46 * t47;
t33 = t43 * t52 + t50;
t53 = t45 * t43;
t29 = atan2(-t33, t53);
t26 = sin(t29);
t27 = cos(t29);
t25 = -t26 * t33 + t27 * t53;
t24 = 0.1e1 / t25 ^ 2;
t51 = t48 * t43;
t35 = -t46 * t44 + t47 * t51;
t58 = t24 * t35;
t56 = t27 * t33;
t55 = t35 ^ 2 * t24;
t38 = 0.1e1 / t43;
t54 = t38 / t45;
t36 = t46 * t43 + t47 * t50;
t32 = 0.1e1 / t36 ^ 2;
t49 = t48 ^ 2 * t59 * t32;
t41 = 0.1e1 / t59;
t31 = 0.1e1 / t36;
t30 = 0.1e1 / (0.1e1 + t49);
t28 = 0.1e1 / (0.1e1 + t33 ^ 2 * t41 / t43 ^ 2);
t23 = 0.1e1 / t25;
t22 = (t33 * t38 * t41 * t47 + t46) * t28;
t21 = 0.1e1 / (0.1e1 + t55);
t1 = [-t35 * t28 * t54, t22, 0, 0, 0, 0; (-t33 * t23 - (-t26 + (t54 * t56 + t26) * t28) * t55) * t21 (t22 * t56 * t58 + (-t48 * t45 * t23 - (t27 * t47 + (-t22 + t46) * t45 * t26) * t58) * t43) * t21, 0, 0, 0, 0; (t46 * t31 + (-t44 * t52 + t51) * t48 * t32) * t45 * t30 (-t31 * t47 * t48 - t44 * t49) * t30, 0, 0, 0, 0;];
Ja_rot  = t1;
