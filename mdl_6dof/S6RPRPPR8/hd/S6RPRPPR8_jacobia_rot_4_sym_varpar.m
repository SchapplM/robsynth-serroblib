% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR8_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobia_rot_4_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:07
% EndTime: 2019-02-26 20:43:07
% DurationCPUTime: 0.05s
% Computational Cost: add. (58->14), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->26)
t27 = sin(qJ(1));
t30 = t27 ^ 2;
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t29 = cos(qJ(1));
t31 = t29 * t28;
t18 = atan2(-t31, t26);
t16 = sin(t18);
t17 = cos(t18);
t14 = -t16 * t31 + t17 * t26;
t13 = 0.1e1 / t14 ^ 2;
t36 = t13 * t28;
t35 = t16 * t26;
t22 = 0.1e1 / t26 ^ 2;
t24 = t28 ^ 2;
t34 = t22 * t24;
t25 = t29 ^ 2;
t33 = 0.1e1 / t30 * t25;
t20 = 0.1e1 / (t25 * t34 + 0.1e1);
t32 = t29 * t20;
t21 = 0.1e1 / t26;
t19 = 0.1e1 / (t22 * t33 + 0.1e1);
t15 = (0.1e1 + t34) * t32;
t12 = 0.1e1 / t14;
t11 = 0.1e1 / (t30 * t24 * t13 + 0.1e1);
t1 = [t27 * t28 * t21 * t20, 0, t15, 0, 0, 0; (-t12 * t31 + (-t17 * t21 * t24 * t32 + (-t20 + 0.1e1) * t28 * t16) * t30 * t36) * t11, 0 (t26 * t12 + (t29 * t35 + t17 * t28 + (-t17 * t31 - t35) * t15) * t36) * t27 * t11, 0, 0, 0; (0.1e1 + t33) * t21 * t19, 0, 0.1e1 / t27 * t22 * t19 * t31, 0, 0, 0;];
Ja_rot  = t1;
