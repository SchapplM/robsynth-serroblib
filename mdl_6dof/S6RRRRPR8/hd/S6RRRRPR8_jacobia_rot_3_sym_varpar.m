% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR8_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_jacobia_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:42
% EndTime: 2019-02-26 22:34:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
t40 = cos(qJ(2));
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t46 = t38 * t37;
t31 = atan2(-t46, -t40);
t29 = sin(t31);
t30 = cos(t31);
t22 = -t29 * t46 - t30 * t40;
t21 = 0.1e1 / t22 ^ 2;
t41 = cos(qJ(1));
t51 = t21 * t41 ^ 2;
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t43 = t41 * t39;
t28 = t38 * t36 + t40 * t43;
t26 = 0.1e1 / t28 ^ 2;
t44 = t41 * t36;
t27 = -t38 * t39 + t40 * t44;
t50 = t26 * t27;
t33 = t37 ^ 2;
t49 = t33 / t40 ^ 2;
t48 = t37 * t41;
t32 = 0.1e1 / (t38 ^ 2 * t49 + 0.1e1);
t47 = t38 * t32;
t45 = t38 * t40;
t42 = t27 ^ 2 * t26 + 0.1e1;
t34 = 0.1e1 / t40;
t25 = 0.1e1 / t28;
t24 = (0.1e1 + t49) * t47;
t23 = 0.1e1 / t42;
t20 = 0.1e1 / t22;
t19 = 0.1e1 / (t33 * t51 + 0.1e1);
t1 = [t34 * t32 * t48, t24, 0, 0, 0, 0; (-t20 * t46 - (-t30 * t33 * t34 * t47 + (t32 - 0.1e1) * t37 * t29) * t37 * t51) * t19 (t40 * t20 - (-t29 * t45 + t30 * t37 + (t29 * t40 - t30 * t46) * t24) * t37 * t21) * t41 * t19, 0, 0, 0, 0; ((-t36 * t45 - t43) * t25 - (-t39 * t45 + t44) * t50) * t23 (-t25 * t36 + t39 * t50) * t23 * t48, t42 * t23, 0, 0, 0;];
Ja_rot  = t1;
