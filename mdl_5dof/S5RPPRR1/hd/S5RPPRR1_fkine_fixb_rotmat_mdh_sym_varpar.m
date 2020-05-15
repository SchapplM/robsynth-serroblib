% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:54
% EndTime: 2019-12-05 17:37:54
% DurationCPUTime: 0.08s
% Computational Cost: add. (56->30), mult. (39->21), div. (0->0), fcn. (67->6), ass. (0->20)
t10 = sin(qJ(4));
t11 = sin(qJ(1));
t23 = t11 * t10;
t13 = cos(qJ(1));
t22 = t13 * t10;
t8 = pkin(5) + 0;
t21 = t11 * pkin(1) + 0;
t20 = pkin(2) + t8;
t3 = t11 * qJ(3);
t19 = t3 + t21;
t18 = t13 * pkin(1) + t11 * qJ(2) + 0;
t17 = pkin(3) + t20;
t16 = t13 * qJ(3) + t18;
t15 = -t13 * qJ(2) + t21;
t14 = -pkin(7) - pkin(6);
t12 = cos(qJ(4));
t9 = qJ(4) + qJ(5);
t2 = cos(t9);
t1 = sin(t9);
t4 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t13, -t11, 0, 0; t11, t13, 0, 0; 0, 0, 1, t8; 0, 0, 0, 1; 0, -t13, t11, t18; 0, -t11, -t13, t15; 1, 0, 0, t8; 0, 0, 0, 1; 0, t11, t13, t16; 0, -t13, t11, t15 + t3; 1, 0, 0, t20; 0, 0, 0, 1; t22, t13 * t12, -t11, -t11 * pkin(6) + t16; t23, t11 * t12, t13, (pkin(6) - qJ(2)) * t13 + t19; t12, -t10, 0, t17; 0, 0, 0, 1; t13 * t1, t13 * t2, -t11, pkin(4) * t22 + t11 * t14 + t16; t11 * t1, t11 * t2, t13, pkin(4) * t23 + (-qJ(2) - t14) * t13 + t19; t2, -t1, 0, t12 * pkin(4) + t17; 0, 0, 0, 1;];
T_ges = t4;
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
