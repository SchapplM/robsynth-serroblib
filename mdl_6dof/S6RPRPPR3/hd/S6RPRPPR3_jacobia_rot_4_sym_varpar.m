% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR3_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobia_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:19
% EndTime: 2019-02-26 20:40:19
% DurationCPUTime: 0.06s
% Computational Cost: add. (153->17), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->27)
t26 = qJ(1) + pkin(9);
t25 = cos(t26);
t32 = t25 ^ 2;
t31 = cos(qJ(3));
t24 = sin(t26);
t30 = sin(qJ(3));
t34 = t24 * t30;
t20 = atan2(-t34, -t31);
t17 = sin(t20);
t18 = cos(t20);
t15 = -t17 * t34 - t18 * t31;
t14 = 0.1e1 / t15 ^ 2;
t38 = t14 * t30;
t37 = t17 * t31;
t22 = t24 ^ 2;
t36 = t22 / t32;
t27 = t30 ^ 2;
t29 = 0.1e1 / t31 ^ 2;
t33 = t27 * t29;
t21 = 0.1e1 / (t22 * t33 + 0.1e1);
t35 = t24 * t21;
t28 = 0.1e1 / t31;
t19 = 0.1e1 / (t29 * t36 + 0.1e1);
t16 = (0.1e1 + t33) * t35;
t13 = 0.1e1 / t15;
t12 = 0.1e1 / (t32 * t27 * t14 + 0.1e1);
t1 = [t25 * t30 * t28 * t21, 0, t16, 0, 0, 0; (-t13 * t34 - (-t18 * t27 * t28 * t35 + (t21 - 0.1e1) * t30 * t17) * t32 * t38) * t12, 0 (t31 * t13 - (-t24 * t37 + t18 * t30 + (-t18 * t34 + t37) * t16) * t38) * t25 * t12, 0, 0, 0; (-0.1e1 - t36) * t28 * t19, 0, -0.1e1 / t25 * t29 * t19 * t34, 0, 0, 0;];
Ja_rot  = t1;
