% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:15
% EndTime: 2020-01-03 11:57:15
% DurationCPUTime: 0.09s
% Computational Cost: add. (128->33), mult. (52->30), div. (0->0), fcn. (85->10), ass. (0->27)
t11 = sin(pkin(9));
t10 = qJ(1) + qJ(2);
t6 = pkin(8) + t10;
t2 = sin(t6);
t29 = t2 * t11;
t3 = cos(t6);
t28 = t3 * t11;
t12 = cos(pkin(9));
t13 = sin(qJ(5));
t27 = t12 * t13;
t15 = cos(qJ(5));
t26 = t12 * t15;
t25 = pkin(5) + 0;
t14 = sin(qJ(1));
t24 = t14 * pkin(1) + 0;
t23 = pkin(6) + t25;
t7 = sin(t10);
t22 = pkin(2) * t7 + t24;
t16 = cos(qJ(1));
t21 = -t16 * pkin(1) + 0;
t5 = qJ(3) + t23;
t20 = pkin(4) * t12 + pkin(7) * t11;
t8 = cos(t10);
t19 = -pkin(2) * t8 + t21;
t18 = t2 * pkin(3) - t3 * qJ(4) + t22;
t17 = -t2 * qJ(4) + t19;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t25; t14, t16, 0, 0; -t16, t14, 0, 0; 0, 0, 0, 1; 0, 0, 1, t23; t7, t8, 0, t24; -t8, t7, 0, t21; 0, 0, 0, 1; 0, 0, 1, t5; t2, t3, 0, t22; -t3, t2, 0, t19; 0, 0, 0, 1; t11, t12, 0, t5; t2 * t12, -t29, -t3, t18; -t3 * t12, t28, -t2, -t3 * pkin(3) + t17; 0, 0, 0, 1; t11 * t15, -t11 * t13, -t12, t11 * pkin(4) - t12 * pkin(7) + t5; -t3 * t13 + t2 * t26, -t3 * t15 - t2 * t27, t29, t20 * t2 + t18; -t2 * t13 - t3 * t26, -t2 * t15 + t3 * t27, -t28, (-pkin(3) - t20) * t3 + t17; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
