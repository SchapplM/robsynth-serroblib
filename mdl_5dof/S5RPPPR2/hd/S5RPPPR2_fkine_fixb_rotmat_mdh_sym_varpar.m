% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR2
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
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:11
% EndTime: 2020-01-03 11:22:11
% DurationCPUTime: 0.18s
% Computational Cost: add. (124->53), mult. (250->61), div. (0->0), fcn. (348->10), ass. (0->37)
t23 = sin(pkin(8));
t24 = sin(pkin(7));
t46 = t24 * t23;
t26 = cos(pkin(8));
t45 = t24 * t26;
t29 = sin(qJ(1));
t44 = t29 * t24;
t27 = cos(pkin(7));
t43 = t29 * t27;
t31 = cos(qJ(1));
t42 = t31 * t24;
t41 = t31 * t27;
t40 = qJ(3) * t24;
t21 = pkin(5) + 0;
t39 = -t29 * qJ(2) + 0;
t38 = t29 * pkin(1) - t31 * qJ(2) + 0;
t37 = t24 * pkin(2) - t27 * qJ(3) + t21;
t36 = pkin(2) * t43 + t29 * t40 + t38;
t35 = pkin(3) * t45 + qJ(4) * t46 + t37;
t10 = -t31 * t23 + t26 * t43;
t9 = t23 * t43 + t31 * t26;
t34 = t10 * pkin(3) + t9 * qJ(4) + t36;
t33 = (-pkin(2) * t27 - pkin(1) - t40) * t31 + t39;
t11 = t23 * t41 - t29 * t26;
t12 = -t29 * t23 - t26 * t41;
t32 = t12 * pkin(3) - t11 * qJ(4) + t33;
t30 = cos(qJ(5));
t28 = sin(qJ(5));
t25 = cos(pkin(9));
t22 = sin(pkin(9));
t8 = -t27 * t22 + t25 * t45;
t7 = t22 * t45 + t27 * t25;
t4 = t12 * t25 - t22 * t42;
t3 = t12 * t22 + t25 * t42;
t2 = t10 * t25 + t22 * t44;
t1 = t10 * t22 - t25 * t44;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t21; t29, t31, 0, 0; -t31, t29, 0, 0; 0, 0, 0, 1; t24, t27, 0, t21; t43, -t44, -t31, t38; -t41, t42, -t29, -t31 * pkin(1) + t39; 0, 0, 0, 1; t45, -t46, -t27, t37; t10, -t9, t44, t36; t12, t11, -t42, t33; 0, 0, 0, 1; t8, -t7, t46, t35; t2, -t1, t9, t34; t4, -t3, -t11, t32; 0, 0, 0, 1; t28 * t46 + t8 * t30, -t8 * t28 + t30 * t46, t7, t8 * pkin(4) + t7 * pkin(6) + t35; t2 * t30 + t9 * t28, -t2 * t28 + t9 * t30, t1, t2 * pkin(4) + t1 * pkin(6) + t34; -t11 * t28 + t4 * t30, -t11 * t30 - t4 * t28, t3, t4 * pkin(4) + t3 * pkin(6) + t32; 0, 0, 0, 1;];
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
