% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP1
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRPRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:48
% EndTime: 2020-01-03 11:58:48
% DurationCPUTime: 0.07s
% Computational Cost: add. (108->31), mult. (33->18), div. (0->0), fcn. (61->8), ass. (0->24)
t14 = sin(qJ(4));
t12 = qJ(1) + qJ(2);
t8 = pkin(8) + t12;
t3 = sin(t8);
t25 = t3 * t14;
t16 = cos(qJ(4));
t4 = cos(t8);
t24 = t4 * t16;
t23 = pkin(5) + 0;
t15 = sin(qJ(1));
t22 = t15 * pkin(1) + 0;
t21 = pkin(6) + t23;
t9 = sin(t12);
t20 = pkin(2) * t9 + t22;
t17 = cos(qJ(1));
t19 = -t17 * pkin(1) + 0;
t7 = qJ(3) + t21;
t10 = cos(t12);
t18 = -pkin(2) * t10 + t19;
t13 = -qJ(5) - pkin(7);
t6 = t16 * pkin(4) + pkin(3);
t2 = t4 * t14;
t1 = t3 * t16;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t23; t15, t17, 0, 0; -t17, t15, 0, 0; 0, 0, 0, 1; 0, 0, 1, t21; t9, t10, 0, t22; -t10, t9, 0, t19; 0, 0, 0, 1; 0, 0, 1, t7; t3, t4, 0, t20; -t4, t3, 0, t18; 0, 0, 0, 1; t14, t16, 0, t7; t1, -t25, -t4, t3 * pkin(3) - t4 * pkin(7) + t20; -t24, t2, -t3, -t4 * pkin(3) - t3 * pkin(7) + t18; 0, 0, 0, 1; t14, t16, 0, t14 * pkin(4) + t7; t1, -t25, -t4, t4 * t13 + t3 * t6 + t20; -t24, t2, -t3, t3 * t13 - t4 * t6 + t18; 0, 0, 0, 1;];
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
