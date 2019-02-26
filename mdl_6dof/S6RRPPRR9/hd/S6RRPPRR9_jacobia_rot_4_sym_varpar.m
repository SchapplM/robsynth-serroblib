% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR9_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:12
% EndTime: 2019-02-26 21:33:12
% DurationCPUTime: 0.07s
% Computational Cost: add. (102->20), mult. (316->53), div. (70->11), fcn. (497->9), ass. (0->33)
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t43 = cos(qJ(2));
t44 = cos(qJ(1));
t47 = cos(pkin(6));
t45 = t44 * t47;
t29 = t41 * t45 + t42 * t43;
t40 = sin(pkin(6));
t48 = t40 * t41;
t27 = atan2(-t29, t48);
t24 = cos(t27);
t52 = t24 * t29;
t46 = t42 * t47;
t31 = -t44 * t41 - t43 * t46;
t34 = 0.1e1 / t40;
t35 = 0.1e1 / t40 ^ 2;
t39 = 0.1e1 / t42 ^ 2;
t51 = 0.1e1 / (t31 ^ 2 * t39 * t35 + 0.1e1) * t34;
t23 = sin(t27);
t22 = -t23 * t29 + t24 * t48;
t21 = 0.1e1 / t22 ^ 2;
t32 = -t41 * t46 + t44 * t43;
t50 = t32 ^ 2 * t21;
t36 = 0.1e1 / t41;
t49 = t34 * t36;
t38 = 0.1e1 / t42;
t37 = 0.1e1 / t41 ^ 2;
t28 = t42 * t41 - t43 * t45;
t25 = 0.1e1 / (t29 ^ 2 * t35 * t37 + 0.1e1);
t20 = 0.1e1 / t22;
t19 = 0.1e1 / (0.1e1 + t50);
t18 = (t29 * t37 * t43 + t28 * t36) * t34 * t25;
t1 = [-t32 * t25 * t49, t18, 0, 0, 0, 0; (-t29 * t20 - (-t23 + (t49 * t52 + t23) * t25) * t50) * t19 (t31 * t20 - (t24 * t40 * t43 + t23 * t28 + (-t23 * t48 - t52) * t18) * t32 * t21) * t19, 0, 0, 0, 0; (-t31 * t39 * t44 + t28 * t38) * t51, -t32 * t38 * t51, 0, 0, 0, 0;];
Ja_rot  = t1;
