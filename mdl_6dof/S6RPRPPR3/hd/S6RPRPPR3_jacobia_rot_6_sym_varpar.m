% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:24
% EndTime: 2019-02-26 20:40:24
% DurationCPUTime: 0.07s
% Computational Cost: add. (196->21), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
t38 = qJ(1) + pkin(9);
t37 = cos(t38);
t56 = t37 ^ 2;
t43 = sin(qJ(3));
t36 = sin(t38);
t45 = cos(qJ(3));
t51 = t36 * t45;
t34 = atan2(-t51, t43);
t32 = sin(t34);
t33 = cos(t34);
t27 = -t32 * t51 + t33 * t43;
t24 = 0.1e1 / t27 ^ 2;
t55 = t24 * t45;
t42 = sin(qJ(6));
t44 = cos(qJ(6));
t47 = t43 * t44;
t31 = -t36 * t42 + t37 * t47;
t29 = 0.1e1 / t31 ^ 2;
t48 = t42 * t43;
t30 = t36 * t44 + t37 * t48;
t54 = t29 * t30;
t53 = t32 * t43;
t41 = t45 ^ 2;
t49 = 0.1e1 / t43 ^ 2 * t41;
t35 = 0.1e1 / (t36 ^ 2 * t49 + 0.1e1);
t52 = t36 * t35;
t50 = t37 * t45;
t46 = t30 ^ 2 * t29 + 0.1e1;
t39 = 0.1e1 / t43;
t28 = 0.1e1 / t31;
t26 = (0.1e1 + t49) * t52;
t25 = 0.1e1 / t46;
t23 = 0.1e1 / t27;
t22 = 0.1e1 / (t56 * t41 * t24 + 0.1e1);
t1 = [-t39 * t35 * t50, 0, t26, 0, 0, 0; (-t23 * t51 - (t33 * t39 * t41 * t52 + (t35 - 0.1e1) * t45 * t32) * t56 * t55) * t22, 0 (-t43 * t23 - (t36 * t53 + t33 * t45 + (-t33 * t51 - t53) * t26) * t55) * t37 * t22, 0, 0, 0; ((-t36 * t48 + t37 * t44) * t28 - (-t36 * t47 - t37 * t42) * t54) * t25, 0 (t28 * t42 - t44 * t54) * t25 * t50, 0, 0, t46 * t25;];
Ja_rot  = t1;
