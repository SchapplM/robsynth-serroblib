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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:19:57
% EndTime: 2019-12-05 18:19:57
% DurationCPUTime: 0.10s
% Computational Cost: add. (132->37), mult. (52->30), div. (0->0), fcn. (85->10), ass. (0->27)
t13 = sin(pkin(9));
t12 = qJ(1) + qJ(2);
t8 = pkin(8) + t12;
t4 = sin(t8);
t31 = t4 * t13;
t5 = cos(t8);
t30 = t5 * t13;
t14 = cos(pkin(9));
t15 = sin(qJ(5));
t29 = t14 * t15;
t17 = cos(qJ(5));
t28 = t14 * t17;
t27 = pkin(5) + 0;
t18 = cos(qJ(1));
t26 = t18 * pkin(1) + 0;
t25 = pkin(6) + t27;
t10 = cos(t12);
t24 = pkin(2) * t10 + t26;
t16 = sin(qJ(1));
t23 = -t16 * pkin(1) + 0;
t7 = qJ(3) + t25;
t22 = pkin(4) * t14 + pkin(7) * t13;
t21 = t5 * pkin(3) + t4 * qJ(4) + t24;
t9 = sin(t12);
t20 = -pkin(2) * t9 + t23;
t19 = t5 * qJ(4) + t20;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t27; -t16, -t18, 0, 0; t18, -t16, 0, 0; 0, 0, 0, 1; 0, 0, 1, t25; -t9, -t10, 0, t23; t10, -t9, 0, t26; 0, 0, 0, 1; 0, 0, 1, t7; -t4, -t5, 0, t20; t5, -t4, 0, t24; 0, 0, 0, 1; t13, t14, 0, t7; -t4 * t14, t31, t5, -t4 * pkin(3) + t19; t5 * t14, -t30, t4, t21; 0, 0, 0, 1; t13 * t17, -t13 * t15, -t14, t13 * pkin(4) - t14 * pkin(7) + t7; t5 * t15 - t4 * t28, t5 * t17 + t4 * t29, -t31, (-pkin(3) - t22) * t4 + t19; t4 * t15 + t5 * t28, t4 * t17 - t5 * t29, t30, t22 * t5 + t21; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
