% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-29 15:51
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 15:50:13
% EndTime: 2019-12-29 15:50:13
% DurationCPUTime: 0.18s
% Computational Cost: add. (83->37), mult. (89->28), div. (0->0), fcn. (141->8), ass. (0->20)
t27 = sin(qJ(1));
t26 = cos(pkin(7));
t25 = sin(pkin(7));
t16 = pkin(5) + 0;
t9 = -qJ(3) + t16;
t20 = cos(qJ(1));
t24 = t20 * pkin(1) + t27 * qJ(2) + 0;
t23 = t20 * pkin(2) + t24;
t22 = t27 * pkin(1) - t20 * qJ(2) + 0;
t21 = t27 * pkin(2) + t22;
t19 = -pkin(6) - qJ(4);
t18 = cos(pkin(8));
t17 = sin(pkin(8));
t15 = pkin(8) + qJ(5);
t8 = cos(t15);
t7 = sin(t15);
t6 = t18 * pkin(4) + pkin(3);
t2 = t20 * t25 - t27 * t26;
t1 = -t20 * t26 - t27 * t25;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t20, -t27, 0, 0; t27, t20, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t20, 0, t27, t24; t27, 0, -t20, t22; 0, 1, 0, t16; 0, 0, 0, 1; -t1, -t2, 0, t23; -t2, t1, 0, t21; 0, 0, -1, t9; 0, 0, 0, 1; -t1 * t18, t1 * t17, t2, -t1 * pkin(3) + t2 * qJ(4) + t23; -t2 * t18, t2 * t17, -t1, -t2 * pkin(3) - t1 * qJ(4) + t21; -t17, -t18, 0, t9; 0, 0, 0, 1; -t1 * t8, t1 * t7, t2, -t1 * t6 - t2 * t19 + t23; -t2 * t8, t2 * t7, -t1, t1 * t19 - t2 * t6 + t21; -t7, -t8, 0, -t17 * pkin(4) + t9; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
