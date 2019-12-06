% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:07
% EndTime: 2019-12-05 17:59:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (73->37), mult. (47->26), div. (0->0), fcn. (79->6), ass. (0->23)
t17 = -pkin(7) - pkin(6);
t13 = sin(qJ(3));
t26 = t13 * pkin(3);
t16 = cos(qJ(1));
t12 = qJ(3) + qJ(4);
t4 = sin(t12);
t25 = t16 * t4;
t5 = cos(t12);
t24 = t16 * t5;
t14 = sin(qJ(1));
t23 = t14 * t13;
t11 = pkin(5) + 0;
t22 = t14 * pkin(1) + 0;
t21 = pkin(2) + t11;
t20 = t16 * pkin(1) + t14 * qJ(2) + 0;
t15 = cos(qJ(3));
t19 = t15 * pkin(3) + t21;
t18 = -t16 * qJ(2) + t22;
t10 = -qJ(5) + t17;
t3 = t14 * t5;
t2 = t14 * t4;
t1 = pkin(4) * t4 + t26;
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t14, 0, 0; t14, t16, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; 0, -t16, t14, t20; 0, -t14, -t16, t18; 1, 0, 0, t11; 0, 0, 0, 1; t23, t14 * t15, t16, t16 * pkin(6) + t20; -t16 * t13, -t16 * t15, t14, t14 * pkin(6) + t18; t15, -t13, 0, t21; 0, 0, 0, 1; t2, t3, t16, pkin(3) * t23 - t16 * t17 + t20; -t25, -t24, t14, -t14 * t17 + (-qJ(2) - t26) * t16 + t22; t5, -t4, 0, t19; 0, 0, 0, 1; t2, t3, t16, t14 * t1 - t16 * t10 + t20; -t25, -t24, t14, -t14 * t10 + (-qJ(2) - t1) * t16 + t22; t5, -t4, 0, pkin(4) * t5 + t19; 0, 0, 0, 1;];
T_ges = t6;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
