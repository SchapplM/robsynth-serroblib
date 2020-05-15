% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:43
% EndTime: 2019-12-05 18:23:44
% DurationCPUTime: 0.09s
% Computational Cost: add. (62->32), mult. (54->31), div. (0->0), fcn. (95->8), ass. (0->25)
t13 = sin(qJ(1));
t9 = qJ(2) + qJ(4);
t7 = sin(t9);
t28 = t13 * t7;
t16 = cos(qJ(1));
t27 = t16 * t7;
t11 = sin(qJ(5));
t26 = t13 * t11;
t12 = sin(qJ(2));
t25 = t13 * t12;
t14 = cos(qJ(5));
t24 = t13 * t14;
t15 = cos(qJ(2));
t4 = t13 * t15;
t23 = t16 * t11;
t22 = t16 * t12;
t21 = t16 * t14;
t5 = t16 * t15;
t17 = pkin(2) + pkin(1);
t20 = t12 * t17 + 0;
t10 = -pkin(3) - qJ(3);
t19 = t16 * t10 + t17 * t4 + 0;
t18 = -t13 * t10 + t17 * t5 + 0;
t8 = cos(t9);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t13, 0, 0; t13, t16, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t5, -t22, t13, 0; t4, -t25, -t16, 0; t12, t15, 0, 0; 0, 0, 0, 1; t5, -t22, t13, pkin(1) * t5 + t13 * qJ(3) + 0; t4, -t25, -t16, pkin(1) * t4 - t16 * qJ(3) + 0; t12, t15, 0, t12 * pkin(1) + 0; 0, 0, 0, 1; t16 * t8, -t27, t13, t18; t13 * t8, -t28, -t16, t19; t7, t8, 0, t20; 0, 0, 0, 1; t8 * t21 + t26, -t8 * t23 + t24, t27, pkin(4) * t27 + t18; t8 * t24 - t23, -t8 * t26 - t21, t28, pkin(4) * t28 + t19; t7 * t14, -t7 * t11, -t8, -t8 * pkin(4) + t20; 0, 0, 0, 1;];
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
