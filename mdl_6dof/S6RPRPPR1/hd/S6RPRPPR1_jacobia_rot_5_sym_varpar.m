% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:11
% EndTime: 2019-02-26 20:39:11
% DurationCPUTime: 0.10s
% Computational Cost: add. (316->22), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->36)
t41 = qJ(3) + pkin(10);
t39 = cos(t41);
t37 = sin(t41);
t42 = qJ(1) + pkin(9);
t38 = sin(t42);
t49 = t38 * t37;
t32 = atan2(-t49, -t39);
t30 = sin(t32);
t31 = cos(t32);
t23 = -t30 * t49 - t31 * t39;
t22 = 0.1e1 / t23 ^ 2;
t40 = cos(t42);
t55 = t22 * t40 ^ 2;
t44 = cos(pkin(11));
t45 = t40 * t44;
t43 = sin(pkin(11));
t48 = t38 * t43;
t29 = t39 * t45 + t48;
t27 = 0.1e1 / t29 ^ 2;
t46 = t40 * t43;
t47 = t38 * t44;
t28 = t39 * t46 - t47;
t54 = t27 * t28;
t53 = t30 * t39;
t34 = t37 ^ 2;
t52 = t34 / t39 ^ 2;
t51 = t37 * t40;
t33 = 0.1e1 / (t38 ^ 2 * t52 + 0.1e1);
t50 = t38 * t33;
t35 = 0.1e1 / t39;
t26 = 0.1e1 / t29;
t25 = 0.1e1 / (t28 ^ 2 * t27 + 0.1e1);
t24 = (0.1e1 + t52) * t50;
t21 = 0.1e1 / t23;
t20 = 0.1e1 / (t34 * t55 + 0.1e1);
t1 = [t35 * t33 * t51, 0, t24, 0, 0, 0; (-t21 * t49 - (-t31 * t34 * t35 * t50 + (t33 - 0.1e1) * t37 * t30) * t37 * t55) * t20, 0 (t39 * t21 - (-t38 * t53 + t31 * t37 + (-t31 * t49 + t53) * t24) * t37 * t22) * t40 * t20, 0, 0, 0; ((-t39 * t48 - t45) * t26 - (-t39 * t47 + t46) * t54) * t25, 0 (-t26 * t43 + t44 * t54) * t25 * t51, 0, 0, 0;];
Ja_rot  = t1;
