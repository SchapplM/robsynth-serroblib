% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR9_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobia_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:12
% EndTime: 2019-02-26 21:33:12
% DurationCPUTime: 0.07s
% Computational Cost: add. (126->20), mult. (316->53), div. (70->11), fcn. (497->9), ass. (0->33)
t43 = sin(qJ(2));
t44 = sin(qJ(1));
t45 = cos(qJ(2));
t46 = cos(qJ(1));
t49 = cos(pkin(6));
t47 = t46 * t49;
t30 = t44 * t43 - t45 * t47;
t42 = sin(pkin(6));
t50 = t42 * t45;
t27 = atan2(-t30, -t50);
t26 = cos(t27);
t54 = t26 * t30;
t48 = t44 * t49;
t34 = -t43 * t48 + t46 * t45;
t36 = 0.1e1 / t42;
t37 = 0.1e1 / t42 ^ 2;
t39 = 0.1e1 / t44 ^ 2;
t53 = 0.1e1 / (t34 ^ 2 * t39 * t37 + 0.1e1) * t36;
t25 = sin(t27);
t24 = -t25 * t30 - t26 * t50;
t23 = 0.1e1 / t24 ^ 2;
t33 = t46 * t43 + t45 * t48;
t52 = t33 ^ 2 * t23;
t40 = 0.1e1 / t45;
t51 = t36 * t40;
t41 = 0.1e1 / t45 ^ 2;
t38 = 0.1e1 / t44;
t32 = t43 * t47 + t44 * t45;
t28 = 0.1e1 / (t30 ^ 2 * t37 * t41 + 0.1e1);
t22 = 0.1e1 / t24;
t21 = 0.1e1 / (0.1e1 + t52);
t20 = (t30 * t41 * t43 + t32 * t40) * t36 * t28;
t1 = [t33 * t28 * t51, t20, 0, 0, 0, 0; (-t30 * t22 - (-t25 + (-t51 * t54 + t25) * t28) * t52) * t21 (t34 * t22 - (t26 * t42 * t43 - t25 * t32 + (t25 * t50 - t54) * t20) * t33 * t23) * t21, 0, 0, 0, 0; (-t34 * t39 * t46 - t32 * t38) * t53, -t33 * t38 * t53, 0, 0, 0, 0;];
Ja_rot  = t1;
