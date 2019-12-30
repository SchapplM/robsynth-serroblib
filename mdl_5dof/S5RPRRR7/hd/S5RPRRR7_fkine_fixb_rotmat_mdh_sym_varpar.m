% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-29 17:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 17:45:03
% EndTime: 2019-12-29 17:45:03
% DurationCPUTime: 0.22s
% Computational Cost: add. (128->42), mult. (83->50), div. (0->0), fcn. (125->10), ass. (0->31)
t16 = sin(qJ(4));
t14 = qJ(1) + pkin(9);
t7 = sin(t14);
t35 = t7 * t16;
t15 = qJ(4) + qJ(5);
t10 = sin(t15);
t20 = cos(qJ(3));
t34 = t10 * t20;
t11 = cos(t15);
t33 = t11 * t20;
t32 = t16 * t20;
t19 = cos(qJ(4));
t31 = t19 * t20;
t30 = pkin(5) + 0;
t18 = sin(qJ(1));
t29 = t18 * pkin(1) + 0;
t21 = cos(qJ(1));
t28 = t21 * pkin(1) + 0;
t27 = t7 * pkin(2) + t29;
t9 = qJ(2) + t30;
t8 = cos(t14);
t26 = t8 * pkin(2) + t7 * pkin(6) + t28;
t17 = sin(qJ(3));
t25 = pkin(3) * t20 + pkin(7) * t17;
t22 = -pkin(8) - pkin(7);
t6 = t19 * pkin(4) + pkin(3);
t24 = -t17 * t22 + t20 * t6;
t23 = -t8 * pkin(6) + t27;
t2 = t8 * t17;
t1 = t7 * t17;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t18, 0, 0; t18, t21, 0, 0; 0, 0, 1, t30; 0, 0, 0, 1; t8, -t7, 0, t28; t7, t8, 0, t29; 0, 0, 1, t9; 0, 0, 0, 1; t8 * t20, -t2, t7, t26; t7 * t20, -t1, -t8, t23; t17, t20, 0, t9; 0, 0, 0, 1; t8 * t31 + t35, t7 * t19 - t8 * t32, t2, t25 * t8 + t26; -t8 * t16 + t7 * t31, -t8 * t19 - t7 * t32, t1, t25 * t7 + t23; t17 * t19, -t17 * t16, -t20, t17 * pkin(3) - t20 * pkin(7) + t9; 0, 0, 0, 1; t7 * t10 + t8 * t33, t7 * t11 - t8 * t34, t2, pkin(4) * t35 + t24 * t8 + t26; -t8 * t10 + t7 * t33, -t8 * t11 - t7 * t34, t1, (-pkin(4) * t16 - pkin(6)) * t8 + t24 * t7 + t27; t17 * t11, -t17 * t10, -t20, t17 * t6 + t20 * t22 + t9; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
