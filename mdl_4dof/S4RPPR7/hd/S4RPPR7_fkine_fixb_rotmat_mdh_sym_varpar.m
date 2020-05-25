% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4RPPR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:32
% EndTime: 2019-12-31 16:41:32
% DurationCPUTime: 0.07s
% Computational Cost: add. (43->26), mult. (29->20), div. (0->0), fcn. (53->6), ass. (0->15)
t11 = sin(qJ(1));
t8 = sin(pkin(6));
t17 = t11 * t8;
t7 = pkin(4) + 0;
t16 = t11 * pkin(1) + 0;
t15 = pkin(2) + t7;
t12 = cos(qJ(1));
t14 = t12 * pkin(1) + t11 * qJ(2) + 0;
t13 = -t12 * qJ(2) + t16;
t10 = -pkin(5) - qJ(3);
t9 = cos(pkin(6));
t6 = pkin(6) + qJ(4);
t2 = cos(t6);
t1 = sin(t6);
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t12, -t11, 0, 0; t11, t12, 0, 0; 0, 0, 1, t7; 0, 0, 0, 1; 0, -t12, t11, t14; 0, -t11, -t12, t13; 1, 0, 0, t7; 0, 0, 0, 1; t17, t11 * t9, t12, t12 * qJ(3) + t14; -t12 * t8, -t12 * t9, t11, t11 * qJ(3) + t13; t9, -t8, 0, t15; 0, 0, 0, 1; t11 * t1, t11 * t2, t12, pkin(3) * t17 - t12 * t10 + t14; -t12 * t1, -t12 * t2, t11, -t11 * t10 + (-pkin(3) * t8 - qJ(2)) * t12 + t16; t2, -t1, 0, t9 * pkin(3) + t15; 0, 0, 0, 1;];
T_ges = t3;
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
