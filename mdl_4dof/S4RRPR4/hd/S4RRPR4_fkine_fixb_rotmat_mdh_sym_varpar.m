% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4RRPR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:22
% EndTime: 2019-12-31 17:02:22
% DurationCPUTime: 0.07s
% Computational Cost: add. (64->25), mult. (25->20), div. (0->0), fcn. (49->8), ass. (0->17)
t18 = pkin(4) + 0;
t14 = sin(qJ(1));
t17 = t14 * pkin(1) + 0;
t15 = cos(qJ(1));
t16 = t15 * pkin(1) + 0;
t6 = pkin(5) + t18;
t13 = -pkin(6) - qJ(3);
t12 = cos(pkin(7));
t11 = sin(pkin(7));
t10 = qJ(1) + qJ(2);
t9 = pkin(7) + qJ(4);
t5 = cos(t10);
t4 = sin(t10);
t3 = cos(t9);
t2 = sin(t9);
t1 = t12 * pkin(3) + pkin(2);
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t15, -t14, 0, 0; t14, t15, 0, 0; 0, 0, 1, t18; 0, 0, 0, 1; t5, -t4, 0, t16; t4, t5, 0, t17; 0, 0, 1, t6; 0, 0, 0, 1; t5 * t12, -t5 * t11, t4, t5 * pkin(2) + t4 * qJ(3) + t16; t4 * t12, -t4 * t11, -t5, t4 * pkin(2) - t5 * qJ(3) + t17; t11, t12, 0, t6; 0, 0, 0, 1; t5 * t3, -t5 * t2, t4, t5 * t1 - t4 * t13 + t16; t4 * t3, -t4 * t2, -t5, t4 * t1 + t5 * t13 + t17; t2, t3, 0, t11 * pkin(3) + t6; 0, 0, 0, 1;];
T_ges = t7;
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
