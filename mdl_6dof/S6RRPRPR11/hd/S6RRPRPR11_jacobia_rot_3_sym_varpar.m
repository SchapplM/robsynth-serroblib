% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR11_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobia_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:46
% EndTime: 2019-02-26 21:43:46
% DurationCPUTime: 0.05s
% Computational Cost: add. (82->16), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->27)
t29 = cos(qJ(1));
t25 = t29 ^ 2;
t28 = cos(qJ(2));
t26 = sin(qJ(2));
t27 = sin(qJ(1));
t32 = t27 * t26;
t18 = atan2(-t32, -t28);
t16 = sin(t18);
t17 = cos(t18);
t14 = -t16 * t32 - t17 * t28;
t13 = 0.1e1 / t14 ^ 2;
t38 = t13 * t26;
t37 = t16 * t28;
t21 = t26 ^ 2;
t31 = t28 ^ 2;
t36 = t21 / t31;
t30 = t27 ^ 2;
t35 = 0.1e1 / t30 * t25;
t34 = t26 * t29;
t19 = 0.1e1 / (t30 * t36 + 0.1e1);
t33 = t27 * t19;
t23 = 0.1e1 / t28;
t20 = 0.1e1 / (t31 * t35 + 0.1e1);
t15 = (0.1e1 + t36) * t33;
t12 = 0.1e1 / t14;
t11 = 0.1e1 / (t25 * t21 * t13 + 0.1e1);
t1 = [t23 * t19 * t34, t15, 0, 0, 0, 0; (-t12 * t32 - (-t17 * t21 * t23 * t33 + (t19 - 0.1e1) * t26 * t16) * t25 * t38) * t11 (t28 * t12 - (-t27 * t37 + t17 * t26 + (-t17 * t32 + t37) * t15) * t38) * t29 * t11, 0, 0, 0, 0; (-0.1e1 - t35) * t28 * t20, -0.1e1 / t27 * t20 * t34, 0, 0, 0, 0;];
Ja_rot  = t1;
