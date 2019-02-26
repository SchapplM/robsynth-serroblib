% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:45
% EndTime: 2019-02-26 20:39:45
% DurationCPUTime: 0.06s
% Computational Cost: add. (263->18), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->29)
t38 = qJ(1) + pkin(9);
t36 = cos(t38);
t32 = t36 ^ 2;
t37 = qJ(3) + pkin(10);
t35 = cos(t37);
t33 = sin(t37);
t34 = sin(t38);
t41 = t34 * t33;
t25 = atan2(-t41, -t35);
t23 = sin(t25);
t24 = cos(t25);
t21 = -t23 * t41 - t24 * t35;
t20 = 0.1e1 / t21 ^ 2;
t47 = t20 * t33;
t46 = t23 * t35;
t28 = t33 ^ 2;
t40 = t35 ^ 2;
t45 = t28 / t40;
t39 = t34 ^ 2;
t44 = 0.1e1 / t39 * t32;
t43 = t33 * t36;
t26 = 0.1e1 / (t39 * t45 + 0.1e1);
t42 = t34 * t26;
t30 = 0.1e1 / t35;
t27 = 0.1e1 / (t40 * t44 + 0.1e1);
t22 = (0.1e1 + t45) * t42;
t19 = 0.1e1 / t21;
t18 = 0.1e1 / (t32 * t28 * t20 + 0.1e1);
t1 = [t30 * t26 * t43, 0, t22, 0, 0, 0; (-t19 * t41 - (-t24 * t28 * t30 * t42 + (t26 - 0.1e1) * t33 * t23) * t32 * t47) * t18, 0 (t35 * t19 - (-t34 * t46 + t24 * t33 + (-t24 * t41 + t46) * t22) * t47) * t36 * t18, 0, 0, 0; (-0.1e1 - t44) * t35 * t27, 0, -0.1e1 / t34 * t27 * t43, 0, 0, 0;];
Ja_rot  = t1;
