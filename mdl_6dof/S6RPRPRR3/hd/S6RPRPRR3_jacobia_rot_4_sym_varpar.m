% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR3_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:19
% EndTime: 2019-02-26 20:50:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (194->21), mult. (197->56), div. (47->9), fcn. (297->9), ass. (0->34)
t30 = qJ(1) + pkin(10);
t29 = cos(t30);
t47 = t29 ^ 2;
t37 = cos(qJ(3));
t28 = sin(t30);
t36 = sin(qJ(3));
t42 = t28 * t36;
t26 = atan2(-t42, -t37);
t24 = sin(t26);
t25 = cos(t26);
t18 = -t24 * t42 - t25 * t37;
t16 = 0.1e1 / t18 ^ 2;
t46 = t16 * t36;
t34 = sin(pkin(11));
t35 = cos(pkin(11));
t38 = t35 * t37;
t23 = t28 * t34 + t29 * t38;
t21 = 0.1e1 / t23 ^ 2;
t39 = t34 * t37;
t22 = -t28 * t35 + t29 * t39;
t45 = t21 * t22;
t44 = t24 * t37;
t31 = t36 ^ 2;
t40 = t31 / t37 ^ 2;
t27 = 0.1e1 / (t28 ^ 2 * t40 + 0.1e1);
t43 = t28 * t27;
t41 = t29 * t36;
t32 = 0.1e1 / t37;
t20 = 0.1e1 / t23;
t19 = (0.1e1 + t40) * t43;
t17 = 0.1e1 / (t22 ^ 2 * t21 + 0.1e1);
t15 = 0.1e1 / t18;
t14 = 0.1e1 / (t47 * t31 * t16 + 0.1e1);
t1 = [t32 * t27 * t41, 0, t19, 0, 0, 0; (-t15 * t42 - (-t25 * t31 * t32 * t43 + (t27 - 0.1e1) * t36 * t24) * t47 * t46) * t14, 0 (t37 * t15 - (-t28 * t44 + t25 * t36 + (-t25 * t42 + t44) * t19) * t46) * t29 * t14, 0, 0, 0; ((-t28 * t39 - t29 * t35) * t20 - (-t28 * t38 + t29 * t34) * t45) * t17, 0 (-t20 * t34 + t35 * t45) * t17 * t41, 0, 0, 0;];
Ja_rot  = t1;
