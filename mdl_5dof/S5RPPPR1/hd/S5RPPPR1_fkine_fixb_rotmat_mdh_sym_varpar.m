% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)
% T_c_stack [(5+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPPPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:04
% EndTime: 2022-01-20 09:12:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (128->42), mult. (83->48), div. (0->0), fcn. (125->10), ass. (0->31)
t16 = sin(pkin(9));
t15 = qJ(1) + pkin(7);
t8 = sin(t15);
t35 = t16 * t8;
t19 = cos(pkin(8));
t34 = t8 * t19;
t10 = cos(t15);
t33 = t10 * t19;
t32 = t16 * t19;
t18 = cos(pkin(9));
t31 = t18 * t19;
t30 = pkin(5) + 0;
t21 = sin(qJ(1));
t29 = t21 * pkin(1) + 0;
t22 = cos(qJ(1));
t28 = t22 * pkin(1) + 0;
t27 = t8 * pkin(2) + t29;
t11 = qJ(2) + t30;
t26 = t10 * pkin(2) + t8 * qJ(3) + t28;
t17 = sin(pkin(8));
t20 = -pkin(6) - qJ(4);
t6 = pkin(4) * t18 + pkin(3);
t25 = -t17 * t20 + t19 * t6;
t24 = pkin(3) * t19 + qJ(4) * t17;
t23 = -qJ(3) * t10 + t27;
t14 = pkin(9) + qJ(5);
t9 = cos(t14);
t7 = sin(t14);
t2 = t10 * t17;
t1 = t8 * t17;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t22, -t21, 0, 0; t21, t22, 0, 0; 0, 0, 1, t30; t10, -t8, 0, t28; t8, t10, 0, t29; 0, 0, 1, t11; t33, -t2, t8, t26; t34, -t1, -t10, t23; t17, t19, 0, t11; t10 * t31 + t35, -t10 * t32 + t18 * t8, t2, t24 * t10 + t26; -t10 * t16 + t8 * t31, -t10 * t18 - t8 * t32, t1, t24 * t8 + t23; t17 * t18, -t17 * t16, -t19, pkin(3) * t17 - qJ(4) * t19 + t11; t9 * t33 + t7 * t8, -t7 * t33 + t8 * t9, t2, pkin(4) * t35 + t25 * t10 + t26; -t10 * t7 + t9 * t34, -t10 * t9 - t7 * t34, t1, t25 * t8 + (-pkin(4) * t16 - qJ(3)) * t10 + t27; t17 * t9, -t17 * t7, -t19, t17 * t6 + t19 * t20 + t11;];
Tc_stack = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
