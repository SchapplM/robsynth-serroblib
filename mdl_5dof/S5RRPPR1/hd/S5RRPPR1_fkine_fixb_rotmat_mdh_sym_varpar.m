% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR1
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:37
% EndTime: 2020-01-03 11:55:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (114->32), mult. (33->22), div. (0->0), fcn. (61->10), ass. (0->23)
t13 = qJ(1) + qJ(2);
t24 = pkin(5) + 0;
t17 = sin(qJ(1));
t23 = t17 * pkin(1) + 0;
t22 = pkin(6) + t24;
t9 = sin(t13);
t21 = pkin(2) * t9 + t23;
t18 = cos(qJ(1));
t20 = -t18 * pkin(1) + 0;
t5 = qJ(3) + t22;
t10 = cos(t13);
t19 = -pkin(2) * t10 + t20;
t16 = -pkin(7) - qJ(4);
t15 = cos(pkin(9));
t14 = sin(pkin(9));
t12 = pkin(9) + qJ(5);
t8 = pkin(8) + t13;
t7 = cos(t12);
t6 = sin(t12);
t3 = t15 * pkin(4) + pkin(3);
t2 = cos(t8);
t1 = sin(t8);
t4 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t24; t17, t18, 0, 0; -t18, t17, 0, 0; 0, 0, 0, 1; 0, 0, 1, t22; t9, t10, 0, t23; -t10, t9, 0, t20; 0, 0, 0, 1; 0, 0, 1, t5; t1, t2, 0, t21; -t2, t1, 0, t19; 0, 0, 0, 1; t14, t15, 0, t5; t1 * t15, -t1 * t14, -t2, t1 * pkin(3) - t2 * qJ(4) + t21; -t2 * t15, t2 * t14, -t1, -t2 * pkin(3) - t1 * qJ(4) + t19; 0, 0, 0, 1; t6, t7, 0, t14 * pkin(4) + t5; t1 * t7, -t1 * t6, -t2, t1 * t3 + t2 * t16 + t21; -t2 * t7, t2 * t6, -t1, t1 * t16 - t2 * t3 + t19; 0, 0, 0, 1;];
T_ges = t4;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
