% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:29
% EndTime: 2019-02-26 20:25:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (316->22), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->36)
t38 = pkin(10) + qJ(4);
t36 = cos(t38);
t34 = sin(t38);
t39 = qJ(1) + pkin(9);
t35 = sin(t39);
t46 = t35 * t34;
t29 = atan2(-t46, -t36);
t27 = sin(t29);
t28 = cos(t29);
t20 = -t27 * t46 - t28 * t36;
t19 = 0.1e1 / t20 ^ 2;
t37 = cos(t39);
t52 = t19 * t37 ^ 2;
t41 = cos(pkin(11));
t42 = t37 * t41;
t40 = sin(pkin(11));
t45 = t35 * t40;
t26 = t36 * t42 + t45;
t24 = 0.1e1 / t26 ^ 2;
t43 = t37 * t40;
t44 = t35 * t41;
t25 = t36 * t43 - t44;
t51 = t24 * t25;
t50 = t27 * t36;
t31 = t34 ^ 2;
t49 = t31 / t36 ^ 2;
t48 = t34 * t37;
t30 = 0.1e1 / (t35 ^ 2 * t49 + 0.1e1);
t47 = t35 * t30;
t32 = 0.1e1 / t36;
t23 = 0.1e1 / t26;
t22 = 0.1e1 / (t25 ^ 2 * t24 + 0.1e1);
t21 = (0.1e1 + t49) * t47;
t18 = 0.1e1 / t20;
t17 = 0.1e1 / (t31 * t52 + 0.1e1);
t1 = [t32 * t30 * t48, 0, 0, t21, 0, 0; (-t18 * t46 - (-t28 * t31 * t32 * t47 + (t30 - 0.1e1) * t34 * t27) * t34 * t52) * t17, 0, 0 (t36 * t18 - (-t35 * t50 + t28 * t34 + (-t28 * t46 + t50) * t21) * t34 * t19) * t37 * t17, 0, 0; ((-t36 * t45 - t42) * t23 - (-t36 * t44 + t43) * t51) * t22, 0, 0 (-t23 * t40 + t41 * t51) * t22 * t48, 0, 0;];
Ja_rot  = t1;
