% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-29 18:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRR13_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 18:01:34
% EndTime: 2019-12-29 18:01:34
% DurationCPUTime: 0.21s
% Computational Cost: add. (79->45), mult. (85->45), div. (0->0), fcn. (127->8), ass. (0->30)
t12 = sin(qJ(4));
t14 = sin(qJ(1));
t33 = t14 * t12;
t13 = sin(qJ(3));
t32 = t14 * t13;
t15 = cos(qJ(4));
t31 = t14 * t15;
t16 = cos(qJ(3));
t30 = t14 * t16;
t17 = cos(qJ(1));
t29 = t17 * t12;
t28 = t17 * t13;
t27 = t17 * t15;
t10 = pkin(5) + 0;
t26 = t14 * pkin(1) + 0;
t25 = pkin(2) + t10;
t6 = t14 * pkin(6);
t24 = t6 + t26;
t23 = t17 * pkin(1) + t14 * qJ(2) + 0;
t22 = t17 * pkin(6) + t23;
t21 = pkin(3) * t13 - pkin(7) * t16;
t18 = -pkin(8) - pkin(7);
t2 = t15 * pkin(4) + pkin(3);
t20 = t13 * t2 + t16 * t18;
t19 = -t17 * qJ(2) + t26;
t11 = qJ(4) + qJ(5);
t4 = cos(t11);
t3 = sin(t11);
t1 = t17 * t16;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t14, 0, 0; t14, t17, 0, 0; 0, 0, 1, t10; 0, 0, 0, 1; 0, -t17, t14, t23; 0, -t14, -t17, t19; 1, 0, 0, t10; 0, 0, 0, 1; t32, t30, t17, t22; -t28, -t1, t14, t19 + t6; t16, -t13, 0, t25; 0, 0, 0, 1; t13 * t31 + t29, -t12 * t32 + t27, -t30, t21 * t14 + t22; -t13 * t27 + t33, t12 * t28 + t31, t1, (-qJ(2) - t21) * t17 + t24; t16 * t15, -t16 * t12, t13, t16 * pkin(3) + t13 * pkin(7) + t25; 0, 0, 0, 1; t17 * t3 + t4 * t32, t17 * t4 - t3 * t32, -t30, pkin(4) * t29 + t20 * t14 + t22; t14 * t3 - t4 * t28, t14 * t4 + t3 * t28, t1, pkin(4) * t33 + (-qJ(2) - t20) * t17 + t24; t16 * t4, -t16 * t3, t13, -t13 * t18 + t16 * t2 + t25; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
