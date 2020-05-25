% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:24
% EndTime: 2019-12-05 16:43:24
% DurationCPUTime: 0.09s
% Computational Cost: add. (104->35), mult. (41->26), div. (0->0), fcn. (73->8), ass. (0->24)
t15 = pkin(8) + qJ(2);
t5 = sin(t15);
t16 = qJ(3) + qJ(4);
t8 = sin(t16);
t27 = t5 * t8;
t6 = cos(t15);
t26 = t6 * t8;
t21 = -pkin(7) - pkin(6);
t20 = cos(qJ(3));
t4 = t20 * pkin(3) + pkin(2);
t17 = sin(pkin(8));
t25 = t17 * pkin(1) + 0;
t18 = cos(pkin(8));
t24 = t18 * pkin(1) + 0;
t23 = qJ(1) + 0;
t7 = pkin(5) + t23;
t19 = sin(qJ(3));
t22 = t19 * pkin(3) + t7;
t14 = -qJ(5) + t21;
t9 = cos(t16);
t3 = t6 * t9;
t2 = t5 * t9;
t1 = pkin(4) * t9 + t4;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t17, 0, 0; t17, t18, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t6, -t5, 0, t24; t5, t6, 0, t25; 0, 0, 1, t7; 0, 0, 0, 1; t6 * t20, -t6 * t19, t5, t6 * pkin(2) + t5 * pkin(6) + t24; t5 * t20, -t5 * t19, -t6, t5 * pkin(2) - t6 * pkin(6) + t25; t19, t20, 0, t7; 0, 0, 0, 1; t3, -t26, t5, -t5 * t21 + t6 * t4 + t24; t2, -t27, -t6, t6 * t21 + t5 * t4 + t25; t8, t9, 0, t22; 0, 0, 0, 1; t3, -t26, t5, t6 * t1 - t5 * t14 + t24; t2, -t27, -t6, t5 * t1 + t6 * t14 + t25; t8, t9, 0, pkin(4) * t8 + t22; 0, 0, 0, 1;];
T_ges = t10;
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
