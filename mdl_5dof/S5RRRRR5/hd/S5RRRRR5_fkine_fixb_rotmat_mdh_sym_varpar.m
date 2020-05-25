% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:12:58
% EndTime: 2020-01-03 12:12:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (114->32), mult. (33->22), div. (0->0), fcn. (61->10), ass. (0->23)
t13 = qJ(1) + qJ(2);
t24 = pkin(5) + 0;
t15 = sin(qJ(1));
t23 = t15 * pkin(1) + 0;
t22 = pkin(6) + t24;
t7 = sin(t13);
t21 = pkin(2) * t7 + t23;
t5 = pkin(7) + t22;
t17 = cos(qJ(1));
t20 = -t17 * pkin(1) + 0;
t9 = cos(t13);
t19 = -pkin(2) * t9 + t20;
t18 = -pkin(9) - pkin(8);
t16 = cos(qJ(4));
t14 = sin(qJ(4));
t12 = qJ(4) + qJ(5);
t10 = qJ(3) + t13;
t8 = cos(t12);
t6 = sin(t12);
t4 = t16 * pkin(4) + pkin(3);
t3 = cos(t10);
t2 = sin(t10);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t24; t15, t17, 0, 0; -t17, t15, 0, 0; 0, 0, 0, 1; 0, 0, 1, t22; t7, t9, 0, t23; -t9, t7, 0, t20; 0, 0, 0, 1; 0, 0, 1, t5; t2, t3, 0, t21; -t3, t2, 0, t19; 0, 0, 0, 1; t14, t16, 0, t5; t2 * t16, -t2 * t14, -t3, t2 * pkin(3) - t3 * pkin(8) + t21; -t3 * t16, t3 * t14, -t2, -t3 * pkin(3) - t2 * pkin(8) + t19; 0, 0, 0, 1; t6, t8, 0, t14 * pkin(4) + t5; t2 * t8, -t2 * t6, -t3, t3 * t18 + t2 * t4 + t21; -t3 * t8, t3 * t6, -t2, t2 * t18 - t3 * t4 + t19; 0, 0, 0, 1;];
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
