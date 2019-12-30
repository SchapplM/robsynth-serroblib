% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-29 20:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 20:57:25
% EndTime: 2019-12-29 20:57:26
% DurationCPUTime: 0.29s
% Computational Cost: add. (120->55), mult. (117->64), div. (0->0), fcn. (168->10), ass. (0->31)
t24 = -pkin(8) - pkin(7);
t18 = sin(qJ(3));
t34 = t18 * pkin(3);
t21 = cos(qJ(3));
t7 = t21 * pkin(3) + pkin(2);
t20 = sin(qJ(1));
t33 = t20 * t18;
t22 = cos(qJ(2));
t32 = t20 * t22;
t23 = cos(qJ(1));
t31 = t23 * t22;
t17 = qJ(3) + qJ(4);
t15 = pkin(5) + 0;
t30 = t20 * pkin(1) + 0;
t29 = t23 * pkin(1) + t20 * pkin(6) + 0;
t19 = sin(qJ(2));
t28 = pkin(2) * t22 + pkin(7) * t19;
t9 = cos(t17);
t1 = pkin(4) * t9 + t7;
t16 = -pkin(9) + t24;
t27 = t1 * t22 - t16 * t19;
t26 = -t19 * t24 + t22 * t7;
t25 = -t23 * pkin(6) + t30;
t10 = qJ(5) + t17;
t8 = sin(t17);
t6 = cos(t10);
t5 = sin(t10);
t4 = t23 * t19;
t3 = t20 * t19;
t2 = pkin(4) * t8 + t34;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t20, 0, 0; t20, t23, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t31, -t4, t20, t29; t32, -t3, -t23, t25; t19, t22, 0, t15; 0, 0, 0, 1; t21 * t31 + t33, -t18 * t31 + t20 * t21, t4, t28 * t23 + t29; -t23 * t18 + t21 * t32, -t18 * t32 - t23 * t21, t3, t28 * t20 + t25; t19 * t21, -t19 * t18, -t22, t19 * pkin(2) - t22 * pkin(7) + t15; 0, 0, 0, 1; t20 * t8 + t9 * t31, t20 * t9 - t8 * t31, t4, pkin(3) * t33 + t26 * t23 + t29; -t23 * t8 + t9 * t32, -t23 * t9 - t8 * t32, t3, (-pkin(6) - t34) * t23 + t26 * t20 + t30; t19 * t9, -t19 * t8, -t22, t19 * t7 + t22 * t24 + t15; 0, 0, 0, 1; t20 * t5 + t6 * t31, t20 * t6 - t5 * t31, t4, t20 * t2 + t27 * t23 + t29; -t23 * t5 + t6 * t32, -t23 * t6 - t5 * t32, t3, (-pkin(6) - t2) * t23 + t27 * t20 + t30; t19 * t6, -t19 * t5, -t22, t19 * t1 + t22 * t16 + t15; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
