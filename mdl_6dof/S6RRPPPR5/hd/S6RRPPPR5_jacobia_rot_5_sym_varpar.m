% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPPR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:13
% EndTime: 2019-02-26 21:24:13
% DurationCPUTime: 0.11s
% Computational Cost: add. (102->20), mult. (329->52), div. (61->11), fcn. (494->9), ass. (0->34)
t46 = sin(qJ(2));
t60 = t46 ^ 2;
t45 = cos(pkin(9));
t44 = sin(pkin(9));
t49 = cos(qJ(1));
t52 = t49 * t44;
t47 = sin(qJ(1));
t48 = cos(qJ(2));
t53 = t47 * t48;
t34 = t45 * t53 - t52;
t54 = t46 * t45;
t30 = atan2(-t34, t54);
t27 = sin(t30);
t28 = cos(t30);
t26 = -t27 * t34 + t28 * t54;
t25 = 0.1e1 / t26 ^ 2;
t51 = t49 * t45;
t37 = t47 * t44 + t48 * t51;
t59 = t25 * t37;
t57 = t28 * t34;
t56 = t37 ^ 2 * t25;
t39 = 0.1e1 / t45;
t55 = t39 / t46;
t36 = -t47 * t45 + t48 * t52;
t33 = 0.1e1 / t36 ^ 2;
t50 = t49 ^ 2 * t60 * t33;
t42 = 0.1e1 / t60;
t32 = 0.1e1 / t36;
t31 = 0.1e1 / (0.1e1 + t50);
t29 = 0.1e1 / (0.1e1 + t34 ^ 2 * t42 / t45 ^ 2);
t24 = 0.1e1 / t26;
t23 = (t34 * t39 * t42 * t48 + t47) * t29;
t22 = 0.1e1 / (0.1e1 + t56);
t1 = [-t37 * t29 * t55, t23, 0, 0, 0, 0; (-t34 * t24 - (-t27 + (t55 * t57 + t27) * t29) * t56) * t22 (t23 * t57 * t59 + (-t49 * t46 * t24 - (t28 * t48 + (-t23 + t47) * t27 * t46) * t59) * t45) * t22, 0, 0, 0, 0; (-t47 * t32 - (-t44 * t53 - t51) * t49 * t33) * t46 * t31 (t32 * t48 * t49 + t44 * t50) * t31, 0, 0, 0, 0;];
Ja_rot  = t1;
