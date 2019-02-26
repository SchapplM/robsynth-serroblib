% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:31:27
% EndTime: 2019-02-26 20:31:27
% DurationCPUTime: 0.10s
% Computational Cost: add. (195->20), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
t36 = qJ(1) + pkin(9);
t34 = sin(t36);
t54 = t34 ^ 2;
t41 = sin(qJ(4));
t35 = cos(t36);
t43 = cos(qJ(4));
t48 = t35 * t43;
t32 = atan2(-t48, t41);
t30 = sin(t32);
t31 = cos(t32);
t24 = -t30 * t48 + t31 * t41;
t22 = 0.1e1 / t24 ^ 2;
t53 = t22 * t43;
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t45 = t41 * t42;
t29 = t34 * t45 + t35 * t40;
t27 = 0.1e1 / t29 ^ 2;
t46 = t40 * t41;
t28 = t34 * t46 - t35 * t42;
t52 = t27 * t28;
t51 = t30 * t41;
t50 = t34 * t43;
t39 = t43 ^ 2;
t47 = 0.1e1 / t41 ^ 2 * t39;
t33 = 0.1e1 / (t35 ^ 2 * t47 + 0.1e1);
t49 = t35 * t33;
t44 = t28 ^ 2 * t27 + 0.1e1;
t37 = 0.1e1 / t41;
t26 = 0.1e1 / t29;
t25 = (0.1e1 + t47) * t49;
t23 = 0.1e1 / t44;
t21 = 0.1e1 / t24;
t20 = 0.1e1 / (t54 * t39 * t22 + 0.1e1);
t1 = [t37 * t33 * t50, 0, 0, t25, 0, 0; (-t21 * t48 + (-t31 * t37 * t39 * t49 + (-t33 + 0.1e1) * t43 * t30) * t54 * t53) * t20, 0, 0 (t41 * t21 + (t35 * t51 + t31 * t43 + (-t31 * t48 - t51) * t25) * t53) * t34 * t20, 0, 0; ((t34 * t42 + t35 * t46) * t26 - (-t34 * t40 + t35 * t45) * t52) * t23, 0, 0 (t26 * t40 - t42 * t52) * t23 * t50, t44 * t23, 0;];
Ja_rot  = t1;
