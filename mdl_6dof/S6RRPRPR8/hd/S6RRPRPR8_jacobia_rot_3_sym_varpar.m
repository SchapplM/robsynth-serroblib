% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR8_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobia_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:57
% EndTime: 2019-02-26 21:41:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (99->20), mult. (197->54), div. (47->9), fcn. (297->9), ass. (0->32)
t34 = cos(qJ(2));
t32 = sin(qJ(2));
t33 = sin(qJ(1));
t39 = t33 * t32;
t25 = atan2(-t39, -t34);
t23 = sin(t25);
t24 = cos(t25);
t16 = -t23 * t39 - t24 * t34;
t15 = 0.1e1 / t16 ^ 2;
t35 = cos(qJ(1));
t44 = t15 * t35 ^ 2;
t30 = sin(pkin(10));
t31 = cos(pkin(10));
t36 = t35 * t31;
t22 = t33 * t30 + t34 * t36;
t20 = 0.1e1 / t22 ^ 2;
t37 = t35 * t30;
t21 = -t33 * t31 + t34 * t37;
t43 = t20 * t21;
t27 = t32 ^ 2;
t42 = t27 / t34 ^ 2;
t41 = t32 * t35;
t26 = 0.1e1 / (t33 ^ 2 * t42 + 0.1e1);
t40 = t33 * t26;
t38 = t33 * t34;
t28 = 0.1e1 / t34;
t19 = 0.1e1 / t22;
t18 = (0.1e1 + t42) * t40;
t17 = 0.1e1 / (t21 ^ 2 * t20 + 0.1e1);
t14 = 0.1e1 / t16;
t13 = 0.1e1 / (t27 * t44 + 0.1e1);
t1 = [t28 * t26 * t41, t18, 0, 0, 0, 0; (-t14 * t39 - (-t24 * t27 * t28 * t40 + (t26 - 0.1e1) * t32 * t23) * t32 * t44) * t13 (t34 * t14 - (-t23 * t38 + t24 * t32 + (t23 * t34 - t24 * t39) * t18) * t32 * t15) * t35 * t13, 0, 0, 0, 0; ((-t30 * t38 - t36) * t19 - (-t31 * t38 + t37) * t43) * t17 (-t19 * t30 + t31 * t43) * t17 * t41, 0, 0, 0, 0;];
Ja_rot  = t1;
