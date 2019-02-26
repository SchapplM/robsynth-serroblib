% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobia_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:35
% EndTime: 2019-02-26 20:28:35
% DurationCPUTime: 0.05s
% Computational Cost: add. (36->14), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->27)
t30 = cos(qJ(1));
t26 = t30 ^ 2;
t27 = sin(qJ(4));
t28 = sin(qJ(1));
t29 = cos(qJ(4));
t34 = t28 * t29;
t21 = atan2(t34, t27);
t17 = sin(t21);
t18 = cos(t21);
t15 = t17 * t34 + t18 * t27;
t14 = 0.1e1 / t15 ^ 2;
t39 = t14 * t29;
t38 = t17 * t27;
t25 = t29 ^ 2;
t31 = t27 ^ 2;
t37 = 0.1e1 / t31 * t25;
t32 = t28 ^ 2;
t36 = 0.1e1 / t32 * t26;
t19 = 0.1e1 / (t32 * t37 + 0.1e1);
t35 = t28 * t19;
t33 = t29 * t30;
t22 = 0.1e1 / t27;
t20 = 0.1e1 / (t31 * t36 + 0.1e1);
t16 = (-0.1e1 - t37) * t35;
t13 = 0.1e1 / t15;
t12 = 0.1e1 / (t26 * t25 * t14 + 0.1e1);
t1 = [t22 * t19 * t33, 0, 0, t16, 0, 0; (t13 * t34 + (t18 * t22 * t25 * t35 + (-t19 + 0.1e1) * t29 * t17) * t26 * t39) * t12, 0, 0 (t27 * t13 + (-t28 * t38 + t18 * t29 + (t18 * t34 - t38) * t16) * t39) * t30 * t12, 0, 0; (0.1e1 + t36) * t27 * t20, 0, 0, -0.1e1 / t28 * t20 * t33, 0, 0;];
Ja_rot  = t1;
