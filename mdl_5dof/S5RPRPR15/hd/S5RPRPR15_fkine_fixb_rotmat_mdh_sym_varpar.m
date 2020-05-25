% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR15
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRPR15_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:31
% EndTime: 2019-12-31 18:36:31
% DurationCPUTime: 0.10s
% Computational Cost: add. (79->45), mult. (85->45), div. (0->0), fcn. (127->8), ass. (0->28)
t12 = sin(pkin(8));
t16 = sin(qJ(1));
t31 = t16 * t12;
t15 = sin(qJ(3));
t30 = t16 * t15;
t17 = cos(qJ(3));
t29 = t16 * t17;
t18 = cos(qJ(1));
t28 = t18 * t12;
t27 = t18 * t15;
t11 = pkin(5) + 0;
t26 = t16 * pkin(1) + 0;
t25 = pkin(2) + t11;
t6 = t16 * pkin(6);
t24 = t6 + t26;
t23 = t18 * pkin(1) + t16 * qJ(2) + 0;
t22 = t18 * pkin(6) + t23;
t13 = cos(pkin(8));
t1 = t13 * pkin(4) + pkin(3);
t14 = -pkin(7) - qJ(4);
t21 = t1 * t15 + t14 * t17;
t20 = pkin(3) * t15 - qJ(4) * t17;
t19 = -t18 * qJ(2) + t26;
t10 = pkin(8) + qJ(5);
t4 = cos(t10);
t3 = sin(t10);
t2 = t18 * t17;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t16, 0, 0; t16, t18, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; 0, -t18, t16, t23; 0, -t16, -t18, t19; 1, 0, 0, t11; 0, 0, 0, 1; t30, t29, t18, t22; -t27, -t2, t16, t19 + t6; t17, -t15, 0, t25; 0, 0, 0, 1; t13 * t30 + t28, -t12 * t30 + t18 * t13, -t29, t20 * t16 + t22; -t13 * t27 + t31, t12 * t27 + t16 * t13, t2, (-qJ(2) - t20) * t18 + t24; t17 * t13, -t17 * t12, t15, t17 * pkin(3) + t15 * qJ(4) + t25; 0, 0, 0, 1; t18 * t3 + t4 * t30, t18 * t4 - t3 * t30, -t29, pkin(4) * t28 + t21 * t16 + t22; t16 * t3 - t4 * t27, t16 * t4 + t3 * t27, t2, pkin(4) * t31 + (-qJ(2) - t21) * t18 + t24; t17 * t4, -t17 * t3, t15, t17 * t1 - t15 * t14 + t25; 0, 0, 0, 1;];
T_ges = t5;
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
