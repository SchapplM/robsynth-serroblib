% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPPR4
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
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR4_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:58
% EndTime: 2019-02-26 20:40:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (221->21), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
t35 = pkin(9) + qJ(3);
t34 = cos(t35);
t33 = sin(t35);
t38 = sin(qJ(1));
t44 = t38 * t33;
t28 = atan2(-t44, -t34);
t26 = sin(t28);
t27 = cos(t28);
t19 = -t26 * t44 - t27 * t34;
t18 = 0.1e1 / t19 ^ 2;
t39 = cos(qJ(1));
t50 = t18 * t39 ^ 2;
t37 = cos(pkin(10));
t40 = t39 * t37;
t36 = sin(pkin(10));
t43 = t38 * t36;
t25 = t34 * t40 + t43;
t23 = 0.1e1 / t25 ^ 2;
t41 = t39 * t36;
t42 = t38 * t37;
t24 = t34 * t41 - t42;
t49 = t23 * t24;
t48 = t26 * t34;
t30 = t33 ^ 2;
t47 = t30 / t34 ^ 2;
t46 = t33 * t39;
t29 = 0.1e1 / (t38 ^ 2 * t47 + 0.1e1);
t45 = t38 * t29;
t31 = 0.1e1 / t34;
t22 = 0.1e1 / t25;
t21 = 0.1e1 / (t24 ^ 2 * t23 + 0.1e1);
t20 = (0.1e1 + t47) * t45;
t17 = 0.1e1 / t19;
t16 = 0.1e1 / (t30 * t50 + 0.1e1);
t1 = [t31 * t29 * t46, 0, t20, 0, 0, 0; (-t17 * t44 - (-t27 * t30 * t31 * t45 + (t29 - 0.1e1) * t33 * t26) * t33 * t50) * t16, 0 (t34 * t17 - (-t38 * t48 + t27 * t33 + (-t27 * t44 + t48) * t20) * t33 * t18) * t39 * t16, 0, 0, 0; ((-t34 * t43 - t40) * t22 - (-t34 * t42 + t41) * t49) * t21, 0 (-t22 * t36 + t37 * t49) * t21 * t46, 0, 0, 0;];
Ja_rot  = t1;
