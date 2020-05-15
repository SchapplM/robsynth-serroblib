% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:43
% EndTime: 2020-01-03 12:01:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (114->32), mult. (33->22), div. (0->0), fcn. (61->10), ass. (0->23)
t13 = qJ(1) + qJ(2);
t24 = pkin(5) + 0;
t15 = sin(qJ(1));
t23 = t15 * pkin(1) + 0;
t22 = pkin(6) + t24;
t8 = sin(t13);
t21 = pkin(2) * t8 + t23;
t17 = cos(qJ(1));
t20 = -t17 * pkin(1) + 0;
t5 = qJ(3) + t22;
t10 = cos(t13);
t19 = -pkin(2) * t10 + t20;
t18 = -pkin(8) - pkin(7);
t16 = cos(qJ(4));
t14 = sin(qJ(4));
t12 = qJ(4) + qJ(5);
t9 = cos(t12);
t7 = sin(t12);
t6 = pkin(9) + t13;
t4 = t16 * pkin(4) + pkin(3);
t2 = cos(t6);
t1 = sin(t6);
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t24; t15, t17, 0, 0; -t17, t15, 0, 0; 0, 0, 0, 1; 0, 0, 1, t22; t8, t10, 0, t23; -t10, t8, 0, t20; 0, 0, 0, 1; 0, 0, 1, t5; t1, t2, 0, t21; -t2, t1, 0, t19; 0, 0, 0, 1; t14, t16, 0, t5; t1 * t16, -t1 * t14, -t2, t1 * pkin(3) - t2 * pkin(7) + t21; -t2 * t16, t2 * t14, -t1, -t2 * pkin(3) - t1 * pkin(7) + t19; 0, 0, 0, 1; t7, t9, 0, t14 * pkin(4) + t5; t1 * t9, -t1 * t7, -t2, t1 * t4 + t2 * t18 + t21; -t2 * t9, t2 * t7, -t1, t1 * t18 - t2 * t4 + t19; 0, 0, 0, 1;];
T_ges = t3;
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
