% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:46:56
% EndTime: 2019-12-05 17:46:57
% DurationCPUTime: 0.09s
% Computational Cost: add. (79->38), mult. (47->30), div. (0->0), fcn. (79->8), ass. (0->22)
t15 = sin(qJ(3));
t25 = t15 * pkin(3);
t16 = sin(qJ(1));
t24 = t16 * t15;
t14 = -qJ(4) - pkin(6);
t13 = pkin(5) + 0;
t12 = qJ(3) + pkin(8);
t23 = t16 * pkin(1) + 0;
t22 = pkin(2) + t13;
t18 = cos(qJ(1));
t21 = t18 * pkin(1) + t16 * qJ(2) + 0;
t17 = cos(qJ(3));
t20 = t17 * pkin(3) + t22;
t19 = -t18 * qJ(2) + t23;
t11 = -pkin(7) + t14;
t6 = qJ(5) + t12;
t5 = cos(t12);
t4 = sin(t12);
t3 = cos(t6);
t2 = sin(t6);
t1 = pkin(4) * t4 + t25;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t16, 0, 0; t16, t18, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; 0, -t18, t16, t21; 0, -t16, -t18, t19; 1, 0, 0, t13; 0, 0, 0, 1; t24, t16 * t17, t18, t18 * pkin(6) + t21; -t18 * t15, -t18 * t17, t16, t16 * pkin(6) + t19; t17, -t15, 0, t22; 0, 0, 0, 1; t16 * t4, t16 * t5, t18, pkin(3) * t24 - t18 * t14 + t21; -t18 * t4, -t18 * t5, t16, -t16 * t14 + (-qJ(2) - t25) * t18 + t23; t5, -t4, 0, t20; 0, 0, 0, 1; t16 * t2, t16 * t3, t18, t16 * t1 - t18 * t11 + t21; -t18 * t2, -t18 * t3, t16, -t16 * t11 + (-qJ(2) - t1) * t18 + t23; t3, -t2, 0, pkin(4) * t5 + t20; 0, 0, 0, 1;];
T_ges = t7;
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
