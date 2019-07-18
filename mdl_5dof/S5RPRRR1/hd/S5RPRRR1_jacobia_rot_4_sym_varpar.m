% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRRR1_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobia_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobia_rot_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:25
% EndTime: 2019-07-18 13:26:25
% DurationCPUTime: 0.11s
% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
t37 = cos(qJ(3));
t34 = sin(qJ(3));
t35 = sin(qJ(1));
t43 = t35 * t34;
t28 = atan2(-t43, -t37);
t26 = sin(t28);
t27 = cos(t28);
t19 = -t26 * t43 - t27 * t37;
t18 = 0.1e1 / t19 ^ 2;
t38 = cos(qJ(1));
t48 = t18 * t38 ^ 2;
t33 = sin(qJ(4));
t36 = cos(qJ(4));
t40 = t38 * t36;
t25 = t35 * t33 + t37 * t40;
t23 = 0.1e1 / t25 ^ 2;
t41 = t38 * t33;
t24 = -t35 * t36 + t37 * t41;
t47 = t23 * t24;
t30 = t34 ^ 2;
t46 = t30 / t37 ^ 2;
t45 = t34 * t38;
t29 = 0.1e1 / (t35 ^ 2 * t46 + 0.1e1);
t44 = t35 * t29;
t42 = t35 * t37;
t39 = t24 ^ 2 * t23 + 0.1e1;
t31 = 0.1e1 / t37;
t22 = 0.1e1 / t25;
t21 = (0.1e1 + t46) * t44;
t20 = 0.1e1 / t39;
t17 = 0.1e1 / t19;
t16 = 0.1e1 / (t30 * t48 + 0.1e1);
t1 = [t31 * t29 * t45, 0, t21, 0, 0; (-t17 * t43 - (-t27 * t30 * t31 * t44 + (t29 - 0.1e1) * t34 * t26) * t34 * t48) * t16, 0, (t37 * t17 - (-t26 * t42 + t27 * t34 + (t26 * t37 - t27 * t43) * t21) * t34 * t18) * t38 * t16, 0, 0; ((-t33 * t42 - t40) * t22 - (-t36 * t42 + t41) * t47) * t20, 0, (-t22 * t33 + t36 * t47) * t20 * t45, t39 * t20, 0;];
Ja_rot  = t1;
