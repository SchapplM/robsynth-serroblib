% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4PRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 12:00:39
% EndTime: 2019-12-29 12:00:39
% DurationCPUTime: 0.14s
% Computational Cost: add. (64->25), mult. (25->20), div. (0->0), fcn. (49->8), ass. (0->17)
t12 = sin(pkin(6));
t18 = t12 * pkin(1) + 0;
t14 = cos(pkin(6));
t17 = t14 * pkin(1) + 0;
t16 = qJ(1) + 0;
t6 = pkin(4) + t16;
t15 = -pkin(5) - qJ(3);
t13 = cos(pkin(7));
t11 = sin(pkin(7));
t10 = pkin(6) + qJ(2);
t9 = pkin(7) + qJ(4);
t5 = cos(t10);
t4 = cos(t9);
t3 = sin(t10);
t2 = sin(t9);
t1 = t13 * pkin(3) + pkin(2);
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t14, -t12, 0, 0; t12, t14, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t5, -t3, 0, t17; t3, t5, 0, t18; 0, 0, 1, t6; 0, 0, 0, 1; t5 * t13, -t5 * t11, t3, t5 * pkin(2) + t3 * qJ(3) + t17; t3 * t13, -t3 * t11, -t5, t3 * pkin(2) - t5 * qJ(3) + t18; t11, t13, 0, t6; 0, 0, 0, 1; t5 * t4, -t5 * t2, t3, t5 * t1 - t3 * t15 + t17; t3 * t4, -t3 * t2, -t5, t3 * t1 + t5 * t15 + t18; t2, t4, 0, t11 * pkin(3) + t6; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
