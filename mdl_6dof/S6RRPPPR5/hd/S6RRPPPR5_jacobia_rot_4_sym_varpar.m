% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function Ja_rot = S6RRPPPR5_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobia_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:13
% EndTime: 2019-02-26 21:24:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (100->20), mult. (313->54), div. (67->11), fcn. (484->9), ass. (0->34)
t41 = sin(pkin(9));
t42 = cos(pkin(9));
t46 = cos(qJ(1));
t47 = t46 * t42;
t44 = sin(qJ(1));
t45 = cos(qJ(2));
t49 = t44 * t45;
t31 = t41 * t49 + t47;
t43 = sin(qJ(2));
t50 = t43 * t41;
t29 = atan2(-t31, t50);
t26 = sin(t29);
t27 = cos(t29);
t25 = -t26 * t31 + t27 * t50;
t24 = 0.1e1 / t25 ^ 2;
t48 = t46 * t41;
t33 = -t44 * t42 + t45 * t48;
t56 = t24 * t33;
t54 = t27 * t31;
t53 = t33 ^ 2 * t24;
t36 = 0.1e1 / t41;
t37 = 0.1e1 / t43;
t52 = t36 * t37;
t38 = 0.1e1 / t43 ^ 2;
t51 = t38 * t45;
t40 = 0.1e1 / t46 ^ 2;
t39 = 0.1e1 / t46;
t34 = t44 * t41 + t45 * t47;
t30 = 0.1e1 / (t34 ^ 2 * t40 * t38 + 0.1e1);
t28 = 0.1e1 / (0.1e1 + t31 ^ 2 * t38 / t41 ^ 2);
t23 = 0.1e1 / t25;
t22 = (t31 * t36 * t51 + t44) * t28;
t21 = 0.1e1 / (0.1e1 + t53);
t1 = [-t33 * t28 * t52, t22, 0, 0, 0, 0; (-t31 * t23 - (-t26 + (t52 * t54 + t26) * t28) * t53) * t21 (t22 * t54 * t56 + (-t46 * t43 * t23 - (t27 * t45 + (-t22 + t44) * t26 * t43) * t56) * t41) * t21, 0, 0, 0, 0; ((-t42 * t49 + t48) * t39 + t44 * t34 * t40) * t37 * t30 (-t34 * t39 * t51 - t42) * t30, 0, 0, 0, 0;];
Ja_rot  = t1;
