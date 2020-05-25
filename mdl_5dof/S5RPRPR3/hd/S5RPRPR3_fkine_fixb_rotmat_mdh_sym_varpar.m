% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:48
% EndTime: 2020-01-03 11:35:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (128->33), mult. (52->30), div. (0->0), fcn. (85->10), ass. (0->27)
t11 = sin(pkin(9));
t10 = qJ(1) + pkin(8);
t8 = qJ(3) + t10;
t3 = sin(t8);
t29 = t3 * t11;
t4 = cos(t8);
t28 = t4 * t11;
t12 = cos(pkin(9));
t13 = sin(qJ(5));
t27 = t12 * t13;
t15 = cos(qJ(5));
t26 = t12 * t15;
t25 = pkin(5) + 0;
t14 = sin(qJ(1));
t24 = t14 * pkin(1) + 0;
t6 = sin(t10);
t23 = pkin(2) * t6 + t24;
t22 = qJ(2) + t25;
t16 = cos(qJ(1));
t21 = -t16 * pkin(1) + 0;
t5 = pkin(6) + t22;
t20 = pkin(4) * t12 + pkin(7) * t11;
t7 = cos(t10);
t19 = -pkin(2) * t7 + t21;
t18 = t3 * pkin(3) - t4 * qJ(4) + t23;
t17 = -t3 * qJ(4) + t19;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t25; t14, t16, 0, 0; -t16, t14, 0, 0; 0, 0, 0, 1; 0, 0, 1, t22; t6, t7, 0, t24; -t7, t6, 0, t21; 0, 0, 0, 1; 0, 0, 1, t5; t3, t4, 0, t23; -t4, t3, 0, t19; 0, 0, 0, 1; t11, t12, 0, t5; t3 * t12, -t29, -t4, t18; -t4 * t12, t28, -t3, -t4 * pkin(3) + t17; 0, 0, 0, 1; t11 * t15, -t11 * t13, -t12, t11 * pkin(4) - t12 * pkin(7) + t5; -t4 * t13 + t3 * t26, -t4 * t15 - t3 * t27, t29, t20 * t3 + t18; -t3 * t13 - t4 * t26, -t3 * t15 + t4 * t27, -t28, (-pkin(3) - t20) * t4 + t17; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
