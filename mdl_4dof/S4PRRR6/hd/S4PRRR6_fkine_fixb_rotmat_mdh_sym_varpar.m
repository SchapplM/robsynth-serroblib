% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4PRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:34
% EndTime: 2019-12-31 16:34:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (62->37), mult. (75->46), div. (0->0), fcn. (113->8), ass. (0->25)
t11 = sin(pkin(7));
t13 = sin(qJ(3));
t27 = t11 * t13;
t16 = cos(qJ(2));
t26 = t11 * t16;
t12 = cos(pkin(7));
t25 = t12 * t16;
t24 = t13 * t16;
t15 = cos(qJ(3));
t23 = t15 * t16;
t22 = t11 * pkin(1) + 0;
t9 = qJ(1) + 0;
t21 = t12 * pkin(1) + t11 * pkin(4) + 0;
t14 = sin(qJ(2));
t20 = pkin(2) * t16 + pkin(5) * t14;
t17 = -pkin(6) - pkin(5);
t3 = t15 * pkin(3) + pkin(2);
t19 = -t14 * t17 + t16 * t3;
t18 = -t12 * pkin(4) + t22;
t10 = qJ(3) + qJ(4);
t5 = cos(t10);
t4 = sin(t10);
t2 = t12 * t14;
t1 = t11 * t14;
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t12, -t11, 0, 0; t11, t12, 0, 0; 0, 0, 1, t9; 0, 0, 0, 1; t25, -t2, t11, t21; t26, -t1, -t12, t18; t14, t16, 0, t9; 0, 0, 0, 1; t12 * t23 + t27, t11 * t15 - t12 * t24, t2, t20 * t12 + t21; t11 * t23 - t12 * t13, -t11 * t24 - t12 * t15, t1, t20 * t11 + t18; t14 * t15, -t14 * t13, -t16, t14 * pkin(2) - t16 * pkin(5) + t9; 0, 0, 0, 1; t11 * t4 + t5 * t25, t11 * t5 - t4 * t25, t2, pkin(3) * t27 + t19 * t12 + t21; -t12 * t4 + t5 * t26, -t12 * t5 - t4 * t26, t1, (-pkin(3) * t13 - pkin(4)) * t12 + t19 * t11 + t22; t14 * t5, -t14 * t4, -t16, t14 * t3 + t16 * t17 + t9; 0, 0, 0, 1;];
T_ges = t6;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
