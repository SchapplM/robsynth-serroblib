% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRPP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:05
% EndTime: 2019-12-05 16:06:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (111->32), mult. (48->26), div. (0->0), fcn. (80->8), ass. (0->25)
t15 = pkin(7) + qJ(2);
t7 = sin(t15);
t16 = qJ(3) + pkin(8);
t8 = sin(t16);
t30 = t7 * t8;
t9 = cos(t15);
t29 = t9 * t8;
t17 = sin(pkin(7));
t28 = t17 * pkin(1) + 0;
t18 = cos(pkin(7));
t27 = t18 * pkin(1) + 0;
t26 = qJ(1) + 0;
t11 = pkin(5) + t26;
t19 = -qJ(4) - pkin(6);
t21 = cos(qJ(3));
t6 = t21 * pkin(3) + pkin(2);
t25 = t9 * t19 + t7 * t6 + t28;
t20 = sin(qJ(3));
t24 = t20 * pkin(3) + t11;
t10 = cos(t16);
t23 = pkin(4) * t10 + qJ(5) * t8;
t22 = -t7 * t19 + t9 * t6 + t27;
t4 = t9 * t10;
t3 = t7 * t10;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t17, 0, 0; t17, t18, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t9, -t7, 0, t27; t7, t9, 0, t28; 0, 0, 1, t11; 0, 0, 0, 1; t9 * t21, -t9 * t20, t7, t9 * pkin(2) + t7 * pkin(6) + t27; t7 * t21, -t7 * t20, -t9, t7 * pkin(2) - t9 * pkin(6) + t28; t20, t21, 0, t11; 0, 0, 0, 1; t4, -t29, t7, t22; t3, -t30, -t9, t25; t8, t10, 0, t24; 0, 0, 0, 1; t4, t7, t29, t23 * t9 + t22; t3, -t9, t30, t23 * t7 + t25; t8, 0, -t10, t8 * pkin(4) - t10 * qJ(5) + t24; 0, 0, 0, 1;];
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
