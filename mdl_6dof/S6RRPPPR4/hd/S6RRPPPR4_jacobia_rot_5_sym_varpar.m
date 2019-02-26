% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPPR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:35
% EndTime: 2019-02-26 21:23:35
% DurationCPUTime: 0.11s
% Computational Cost: add. (102->19), mult. (329->51), div. (61->11), fcn. (494->9), ass. (0->35)
t47 = cos(qJ(2));
t60 = t47 ^ 2;
t45 = sin(qJ(2));
t43 = sin(pkin(9));
t48 = cos(qJ(1));
t51 = t48 * t43;
t44 = cos(pkin(9));
t46 = sin(qJ(1));
t53 = t46 * t44;
t36 = t45 * t53 + t51;
t52 = t47 * t44;
t30 = atan2(t36, t52);
t27 = sin(t30);
t28 = cos(t30);
t26 = t27 * t36 + t28 * t52;
t25 = 0.1e1 / t26 ^ 2;
t50 = t48 * t44;
t54 = t46 * t43;
t34 = -t45 * t50 + t54;
t59 = t25 * t34;
t57 = t28 * t36;
t56 = t34 ^ 2 * t25;
t38 = 0.1e1 / t44;
t55 = t38 / t47;
t35 = t45 * t51 + t53;
t33 = 0.1e1 / t35 ^ 2;
t49 = t48 ^ 2 * t60 * t33;
t41 = 0.1e1 / t60;
t32 = 0.1e1 / t35;
t31 = 0.1e1 / (0.1e1 + t49);
t29 = 0.1e1 / (0.1e1 + t36 ^ 2 * t41 / t44 ^ 2);
t24 = 0.1e1 / t26;
t23 = (t36 * t38 * t41 * t45 + t46) * t29;
t22 = 0.1e1 / (0.1e1 + t56);
t1 = [-t34 * t29 * t55, t23, 0, 0, 0, 0; (t36 * t24 - (-t27 + (-t55 * t57 + t27) * t29) * t56) * t22 (-t23 * t57 * t59 + (-t48 * t47 * t24 - (-t28 * t45 + (-t23 + t46) * t27 * t47) * t59) * t44) * t22, 0, 0, 0, 0; (t46 * t32 + (-t45 * t54 + t50) * t48 * t33) * t47 * t31 (t32 * t45 * t48 + t43 * t49) * t31, 0, 0, 0, 0;];
Ja_rot  = t1;
