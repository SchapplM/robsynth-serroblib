% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRR4
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4PRRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:26
% EndTime: 2019-12-31 16:32:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (64->25), mult. (25->20), div. (0->0), fcn. (49->8), ass. (0->17)
t11 = sin(pkin(7));
t18 = t11 * pkin(1) + 0;
t12 = cos(pkin(7));
t17 = t12 * pkin(1) + 0;
t16 = qJ(1) + 0;
t4 = pkin(4) + t16;
t15 = -pkin(6) - pkin(5);
t14 = cos(qJ(3));
t13 = sin(qJ(3));
t10 = qJ(3) + qJ(4);
t9 = pkin(7) + qJ(2);
t6 = cos(t10);
t5 = sin(t10);
t3 = cos(t9);
t2 = sin(t9);
t1 = t14 * pkin(3) + pkin(2);
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t12, -t11, 0, 0; t11, t12, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t3, -t2, 0, t17; t2, t3, 0, t18; 0, 0, 1, t4; 0, 0, 0, 1; t3 * t14, -t3 * t13, t2, t3 * pkin(2) + t2 * pkin(5) + t17; t2 * t14, -t2 * t13, -t3, t2 * pkin(2) - t3 * pkin(5) + t18; t13, t14, 0, t4; 0, 0, 0, 1; t3 * t6, -t3 * t5, t2, t3 * t1 - t2 * t15 + t17; t2 * t6, -t2 * t5, -t3, t2 * t1 + t3 * t15 + t18; t5, t6, 0, t13 * pkin(3) + t4; 0, 0, 0, 1;];
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
