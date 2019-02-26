% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR3_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:11
% EndTime: 2019-02-26 21:02:11
% DurationCPUTime: 0.10s
% Computational Cost: add. (218->21), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
t36 = qJ(1) + pkin(10);
t35 = cos(t36);
t54 = t35 ^ 2;
t43 = cos(qJ(3));
t34 = sin(t36);
t41 = sin(qJ(3));
t49 = t34 * t41;
t32 = atan2(-t49, -t43);
t30 = sin(t32);
t31 = cos(t32);
t23 = -t30 * t49 - t31 * t43;
t22 = 0.1e1 / t23 ^ 2;
t53 = t22 * t41;
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t45 = t42 * t43;
t29 = t34 * t40 + t35 * t45;
t27 = 0.1e1 / t29 ^ 2;
t46 = t40 * t43;
t28 = -t34 * t42 + t35 * t46;
t52 = t27 * t28;
t51 = t30 * t43;
t37 = t41 ^ 2;
t47 = t37 / t43 ^ 2;
t33 = 0.1e1 / (t34 ^ 2 * t47 + 0.1e1);
t50 = t34 * t33;
t48 = t35 * t41;
t44 = t28 ^ 2 * t27 + 0.1e1;
t38 = 0.1e1 / t43;
t26 = 0.1e1 / t29;
t25 = (0.1e1 + t47) * t50;
t24 = 0.1e1 / t44;
t21 = 0.1e1 / t23;
t20 = 0.1e1 / (t54 * t37 * t22 + 0.1e1);
t1 = [t38 * t33 * t48, 0, t25, 0, 0, 0; (-t21 * t49 - (-t31 * t37 * t38 * t50 + (t33 - 0.1e1) * t41 * t30) * t54 * t53) * t20, 0 (t43 * t21 - (-t34 * t51 + t31 * t41 + (-t31 * t49 + t51) * t25) * t53) * t35 * t20, 0, 0, 0; ((-t34 * t46 - t35 * t42) * t26 - (-t34 * t45 + t35 * t40) * t52) * t24, 0 (-t26 * t40 + t42 * t52) * t24 * t48, t44 * t24, 0, 0;];
Ja_rot  = t1;
