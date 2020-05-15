% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRPRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:49
% EndTime: 2019-12-31 19:50:49
% DurationCPUTime: 0.08s
% Computational Cost: add. (111->32), mult. (48->26), div. (0->0), fcn. (80->8), ass. (0->25)
t15 = pkin(8) + qJ(4);
t7 = sin(t15);
t16 = qJ(1) + qJ(2);
t9 = sin(t16);
t30 = t9 * t7;
t10 = cos(t16);
t29 = t10 * t7;
t28 = pkin(5) + 0;
t20 = sin(qJ(1));
t27 = t20 * pkin(1) + 0;
t21 = cos(qJ(1));
t26 = t21 * pkin(1) + 0;
t12 = pkin(6) + t28;
t19 = -pkin(7) - qJ(3);
t18 = cos(pkin(8));
t6 = t18 * pkin(3) + pkin(2);
t25 = t10 * t19 + t9 * t6 + t27;
t17 = sin(pkin(8));
t24 = t17 * pkin(3) + t12;
t8 = cos(t15);
t23 = pkin(4) * t8 + qJ(5) * t7;
t22 = t10 * t6 - t9 * t19 + t26;
t4 = t10 * t8;
t3 = t9 * t8;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t20, 0, 0; t20, t21, 0, 0; 0, 0, 1, t28; 0, 0, 0, 1; t10, -t9, 0, t26; t9, t10, 0, t27; 0, 0, 1, t12; 0, 0, 0, 1; t10 * t18, -t10 * t17, t9, t10 * pkin(2) + t9 * qJ(3) + t26; t9 * t18, -t9 * t17, -t10, t9 * pkin(2) - t10 * qJ(3) + t27; t17, t18, 0, t12; 0, 0, 0, 1; t4, -t29, t9, t22; t3, -t30, -t10, t25; t7, t8, 0, t24; 0, 0, 0, 1; t4, t9, t29, t23 * t10 + t22; t3, -t10, t30, t23 * t9 + t25; t7, 0, -t8, t7 * pkin(4) - t8 * qJ(5) + t24; 0, 0, 0, 1;];
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
