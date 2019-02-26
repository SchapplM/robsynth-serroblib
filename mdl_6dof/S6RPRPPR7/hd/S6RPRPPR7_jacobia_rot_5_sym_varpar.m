% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:42:33
% EndTime: 2019-02-26 20:42:33
% DurationCPUTime: 0.06s
% Computational Cost: add. (168->15), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->28)
t35 = sin(qJ(1));
t33 = t35 ^ 2;
t32 = qJ(3) + pkin(9);
t30 = sin(t32);
t31 = cos(t32);
t36 = cos(qJ(1));
t39 = t36 * t31;
t24 = atan2(-t39, t30);
t22 = sin(t24);
t23 = cos(t24);
t20 = -t22 * t39 + t23 * t30;
t19 = 0.1e1 / t20 ^ 2;
t45 = t19 * t31;
t44 = t22 * t30;
t29 = t31 ^ 2;
t37 = t30 ^ 2;
t43 = 0.1e1 / t37 * t29;
t42 = t31 * t35;
t38 = t36 ^ 2;
t41 = t33 / t38;
t25 = 0.1e1 / (t38 * t43 + 0.1e1);
t40 = t36 * t25;
t27 = 0.1e1 / t30;
t26 = 0.1e1 / (t37 * t41 + 0.1e1);
t21 = (0.1e1 + t43) * t40;
t18 = 0.1e1 / t20;
t17 = 0.1e1 / (t33 * t29 * t19 + 0.1e1);
t1 = [t27 * t25 * t42, 0, t21, 0, 0, 0; (-t18 * t39 + (-t23 * t27 * t29 * t40 + (-t25 + 0.1e1) * t31 * t22) * t33 * t45) * t17, 0 (t30 * t18 + (t36 * t44 + t23 * t31 + (-t23 * t39 - t44) * t21) * t45) * t35 * t17, 0, 0, 0; (0.1e1 + t41) * t30 * t26, 0, 0.1e1 / t36 * t26 * t42, 0, 0, 0;];
Ja_rot  = t1;
