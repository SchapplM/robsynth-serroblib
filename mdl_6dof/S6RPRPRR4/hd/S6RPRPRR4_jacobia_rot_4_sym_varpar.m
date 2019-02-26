% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR4_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:44
% EndTime: 2019-02-26 20:50:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (153->17), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->28)
t26 = qJ(1) + pkin(10);
t25 = cos(t26);
t23 = t25 ^ 2;
t31 = cos(qJ(3));
t24 = sin(t26);
t30 = sin(qJ(3));
t36 = t24 * t30;
t20 = atan2(-t36, -t31);
t17 = sin(t20);
t18 = cos(t20);
t15 = -t17 * t36 - t18 * t31;
t14 = 0.1e1 / t15 ^ 2;
t40 = t14 * t30;
t39 = t17 * t31;
t32 = t24 ^ 2;
t38 = 0.1e1 / t32 * t23;
t27 = t30 ^ 2;
t33 = t31 ^ 2;
t34 = t27 / t33;
t21 = 0.1e1 / (t32 * t34 + 0.1e1);
t37 = t24 * t21;
t35 = t25 * t30;
t28 = 0.1e1 / t31;
t19 = 0.1e1 / (t33 * t38 + 0.1e1);
t16 = (0.1e1 + t34) * t37;
t13 = 0.1e1 / t15;
t12 = 0.1e1 / (t23 * t27 * t14 + 0.1e1);
t1 = [t28 * t21 * t35, 0, t16, 0, 0, 0; (-t13 * t36 - (-t18 * t27 * t28 * t37 + (t21 - 0.1e1) * t30 * t17) * t23 * t40) * t12, 0 (t31 * t13 - (-t24 * t39 + t18 * t30 + (-t18 * t36 + t39) * t16) * t40) * t25 * t12, 0, 0, 0; (-0.1e1 - t38) * t31 * t19, 0, -0.1e1 / t24 * t19 * t35, 0, 0, 0;];
Ja_rot  = t1;
