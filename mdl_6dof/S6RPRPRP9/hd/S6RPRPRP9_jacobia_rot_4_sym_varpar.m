% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP9_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobia_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:01
% EndTime: 2019-02-26 20:48:01
% DurationCPUTime: 0.07s
% Computational Cost: add. (76->19), mult. (197->54), div. (47->9), fcn. (297->9), ass. (0->35)
t33 = sin(qJ(1));
t47 = t33 ^ 2;
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t35 = cos(qJ(1));
t36 = t35 * t34;
t25 = atan2(-t36, t32);
t23 = sin(t25);
t24 = cos(t25);
t16 = -t23 * t36 + t24 * t32;
t15 = 0.1e1 / t16 ^ 2;
t46 = t15 * t34;
t30 = sin(pkin(9));
t38 = t35 * t30;
t31 = cos(pkin(9));
t41 = t33 * t31;
t22 = t32 * t41 + t38;
t20 = 0.1e1 / t22 ^ 2;
t37 = t35 * t31;
t42 = t33 * t30;
t21 = t32 * t42 - t37;
t45 = t20 * t21;
t44 = t23 * t32;
t29 = t34 ^ 2;
t43 = 0.1e1 / t32 ^ 2 * t29;
t40 = t33 * t34;
t26 = 0.1e1 / (t35 ^ 2 * t43 + 0.1e1);
t39 = t35 * t26;
t27 = 0.1e1 / t32;
t19 = 0.1e1 / t22;
t18 = (0.1e1 + t43) * t39;
t17 = 0.1e1 / (t21 ^ 2 * t20 + 0.1e1);
t14 = 0.1e1 / t16;
t13 = 0.1e1 / (t47 * t29 * t15 + 0.1e1);
t1 = [t27 * t26 * t40, 0, t18, 0, 0, 0; (-t14 * t36 + (-t24 * t27 * t29 * t39 + (-t26 + 0.1e1) * t34 * t23) * t47 * t46) * t13, 0 (t32 * t14 + (t35 * t44 + t24 * t34 + (-t24 * t36 - t44) * t18) * t46) * t33 * t13, 0, 0, 0; ((t32 * t38 + t41) * t19 - (t32 * t37 - t42) * t45) * t17, 0 (t19 * t30 - t31 * t45) * t17 * t40, 0, 0, 0;];
Ja_rot  = t1;
