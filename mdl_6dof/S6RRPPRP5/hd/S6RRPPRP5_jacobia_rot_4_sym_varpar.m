% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP5_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_jacobia_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:27:36
% EndTime: 2019-02-26 21:27:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (77->20), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->34)
t32 = sin(qJ(2));
t33 = sin(qJ(1));
t34 = cos(qJ(2));
t39 = t33 * t34;
t25 = atan2(-t39, t32);
t23 = sin(t25);
t24 = cos(t25);
t16 = -t23 * t39 + t24 * t32;
t15 = 0.1e1 / t16 ^ 2;
t35 = cos(qJ(1));
t46 = t15 * t35 ^ 2;
t30 = sin(pkin(9));
t37 = t35 * t30;
t31 = cos(pkin(9));
t40 = t33 * t31;
t22 = t32 * t37 + t40;
t20 = 0.1e1 / t22 ^ 2;
t36 = t35 * t31;
t41 = t33 * t30;
t21 = -t32 * t36 + t41;
t45 = t20 * t21;
t44 = t23 * t32;
t29 = t34 ^ 2;
t43 = 0.1e1 / t32 ^ 2 * t29;
t26 = 0.1e1 / (t33 ^ 2 * t43 + 0.1e1);
t42 = t33 * t26;
t38 = t34 * t35;
t27 = 0.1e1 / t32;
t19 = 0.1e1 / t22;
t18 = (0.1e1 + t43) * t42;
t17 = 0.1e1 / (t21 ^ 2 * t20 + 0.1e1);
t14 = 0.1e1 / t16;
t13 = 0.1e1 / (t29 * t46 + 0.1e1);
t1 = [-t27 * t26 * t38, t18, 0, 0, 0, 0; (-t14 * t39 - (t24 * t27 * t29 * t42 + (t26 - 0.1e1) * t34 * t23) * t34 * t46) * t13 (-t32 * t14 - (t33 * t44 + t24 * t34 + (-t24 * t39 - t44) * t18) * t34 * t15) * t35 * t13, 0, 0, 0, 0; ((t32 * t40 + t37) * t19 - (-t32 * t41 + t36) * t45) * t17 (-t19 * t31 - t30 * t45) * t17 * t38, 0, 0, 0, 0;];
Ja_rot  = t1;
