% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR6_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:12
% EndTime: 2019-02-26 19:49:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (89->17), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
t41 = sin(pkin(10));
t42 = sin(pkin(6));
t53 = t41 * t42;
t46 = sin(qJ(2));
t52 = t42 * t46;
t44 = cos(pkin(6));
t51 = t44 * t46;
t48 = cos(qJ(2));
t50 = t44 * t48;
t43 = cos(pkin(10));
t37 = t41 * t50 + t43 * t46;
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t29 = t37 * t45 + t47 * t53;
t27 = 0.1e1 / t29 ^ 2;
t28 = -t37 * t47 + t45 * t53;
t49 = t28 ^ 2 * t27 + 0.1e1;
t40 = 0.1e1 / t46 ^ 2;
t38 = -t41 * t51 + t43 * t48;
t35 = t41 * t48 + t43 * t51;
t34 = t41 * t46 - t43 * t50;
t33 = atan2(-t35, t52);
t31 = cos(t33);
t30 = sin(t33);
t26 = 0.1e1 / t49;
t25 = -t30 * t35 + t31 * t52;
t24 = 0.1e1 / t25 ^ 2;
t22 = (t34 / t46 + t48 * t35 * t40) / t42 / (0.1e1 + t35 ^ 2 / t42 ^ 2 * t40);
t1 = [0, t22, 0, 0, 0, 0; 0 (-t37 / t25 - (t31 * t42 * t48 + t30 * t34 + (-t30 * t52 - t31 * t35) * t22) * t38 * t24) / (t38 ^ 2 * t24 + 0.1e1) 0, 0, 0, 0; 0 (-t47 / t29 - t45 * t28 * t27) * t38 * t26, 0, t49 * t26, 0, 0;];
Ja_rot  = t1;
