% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR5
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:23
% EndTime: 2020-01-03 12:03:23
% DurationCPUTime: 0.12s
% Computational Cost: add. (113->39), mult. (41->30), div. (0->0), fcn. (73->10), ass. (0->23)
t18 = cos(pkin(9));
t4 = t18 * pkin(3) + pkin(2);
t19 = -pkin(7) - qJ(3);
t25 = pkin(5) + 0;
t15 = pkin(9) + qJ(4);
t20 = sin(qJ(1));
t24 = t20 * pkin(1) + 0;
t12 = pkin(6) + t25;
t21 = cos(qJ(1));
t23 = -t21 * pkin(1) + 0;
t17 = sin(pkin(9));
t22 = t17 * pkin(3) + t12;
t16 = qJ(1) + qJ(2);
t14 = -pkin(8) + t19;
t9 = cos(t16);
t8 = sin(t16);
t7 = qJ(5) + t15;
t6 = cos(t15);
t5 = sin(t15);
t3 = cos(t7);
t2 = sin(t7);
t1 = pkin(4) * t6 + t4;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t25; t20, t21, 0, 0; -t21, t20, 0, 0; 0, 0, 0, 1; 0, 0, 1, t12; t8, t9, 0, t24; -t9, t8, 0, t23; 0, 0, 0, 1; t17, t18, 0, t12; t8 * t18, -t8 * t17, -t9, t8 * pkin(2) - t9 * qJ(3) + t24; -t9 * t18, t9 * t17, -t8, -t9 * pkin(2) - t8 * qJ(3) + t23; 0, 0, 0, 1; t5, t6, 0, t22; t8 * t6, -t8 * t5, -t9, t9 * t19 + t8 * t4 + t24; -t9 * t6, t9 * t5, -t8, t8 * t19 - t9 * t4 + t23; 0, 0, 0, 1; t2, t3, 0, pkin(4) * t5 + t22; t8 * t3, -t8 * t2, -t9, t8 * t1 + t9 * t14 + t24; -t9 * t3, t9 * t2, -t8, -t9 * t1 + t8 * t14 + t23; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
