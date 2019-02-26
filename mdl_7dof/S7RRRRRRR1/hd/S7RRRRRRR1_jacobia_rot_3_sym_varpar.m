% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_rot [3x7]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S7RRRRRRR1_jacobia_rot_3_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_rot_3_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_rot_3_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:53
% DurationCPUTime: 0.10s
% Computational Cost: add. (63->18), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
t38 = cos(qJ(2));
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t44 = t36 * t35;
t30 = atan2(t44, t38);
t27 = sin(t30);
t28 = cos(t30);
t20 = t27 * t44 + t28 * t38;
t19 = 0.1e1 / t20 ^ 2;
t39 = cos(qJ(1));
t49 = t19 * t39 ^ 2;
t34 = sin(qJ(3));
t37 = cos(qJ(3));
t41 = t39 * t37;
t26 = -t36 * t34 + t38 * t41;
t24 = 0.1e1 / t26 ^ 2;
t42 = t39 * t34;
t25 = t36 * t37 + t38 * t42;
t48 = t24 * t25;
t31 = t35 ^ 2;
t47 = t31 / t38 ^ 2;
t46 = t35 * t39;
t29 = 0.1e1 / (t36 ^ 2 * t47 + 0.1e1);
t45 = t36 * t29;
t43 = t36 * t38;
t40 = t25 ^ 2 * t24 + 0.1e1;
t32 = 0.1e1 / t38;
t23 = 0.1e1 / t26;
t22 = (0.1e1 + t47) * t45;
t21 = 0.1e1 / t40;
t18 = 0.1e1 / t20;
t17 = 0.1e1 / (t31 * t49 + 0.1e1);
t1 = [t32 * t29 * t46, t22, 0, 0, 0, 0, 0; (t18 * t44 + (t28 * t31 * t32 * t45 + (-t29 + 0.1e1) * t35 * t27) * t35 * t49) * t17 (-t38 * t18 + (t27 * t43 - t28 * t35 + (-t27 * t38 + t28 * t44) * t22) * t35 * t19) * t39 * t17, 0, 0, 0, 0, 0; ((-t34 * t43 + t41) * t23 - (-t37 * t43 - t42) * t48) * t21 (-t23 * t34 + t37 * t48) * t21 * t46, t40 * t21, 0, 0, 0, 0;];
Ja_rot  = t1;
