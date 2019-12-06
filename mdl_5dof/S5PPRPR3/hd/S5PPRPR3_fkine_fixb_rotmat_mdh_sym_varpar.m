% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PPRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:37
% EndTime: 2019-12-05 15:04:38
% DurationCPUTime: 0.13s
% Computational Cost: add. (128->52), mult. (153->63), div. (0->0), fcn. (216->10), ass. (0->36)
t20 = qJ(3) + pkin(9);
t14 = sin(t20);
t21 = sin(pkin(8));
t44 = t21 * t14;
t15 = cos(t20);
t43 = t21 * t15;
t22 = sin(pkin(7));
t10 = t22 * t21;
t23 = cos(pkin(8));
t42 = t22 * t23;
t27 = sin(qJ(3));
t41 = t22 * t27;
t29 = cos(qJ(3));
t40 = t22 * t29;
t24 = cos(pkin(7));
t11 = t24 * t21;
t39 = t24 * t23;
t38 = t24 * t27;
t37 = t24 * t29;
t36 = t22 * pkin(1) + 0;
t19 = qJ(1) + 0;
t35 = t24 * pkin(1) + t22 * qJ(2) + 0;
t13 = t29 * pkin(3) + pkin(2);
t25 = -qJ(4) - pkin(5);
t34 = t21 * t13 + t23 * t25 + t19;
t33 = pkin(2) * t23 + pkin(5) * t21;
t32 = -t24 * qJ(2) + t36;
t31 = pkin(3) * t41 - t25 * t11 + t13 * t39 + t35;
t30 = -t25 * t10 + t13 * t42 + (-pkin(3) * t27 - qJ(2)) * t24 + t36;
t28 = cos(qJ(5));
t26 = sin(qJ(5));
t4 = t22 * t14 + t15 * t39;
t3 = t14 * t39 - t22 * t15;
t2 = -t24 * t14 + t15 * t42;
t1 = t14 * t42 + t24 * t15;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t22, 0, 0; t22, t24, 0, 0; 0, 0, 1, t19; 0, 0, 0, 1; t39, -t11, t22, t35; t42, -t10, -t24, t32; t21, t23, 0, t19; 0, 0, 0, 1; t23 * t37 + t41, -t23 * t38 + t40, t11, t33 * t24 + t35; t23 * t40 - t38, -t23 * t41 - t37, t10, t33 * t22 + t32; t21 * t29, -t21 * t27, -t23, t21 * pkin(2) - t23 * pkin(5) + t19; 0, 0, 0, 1; t4, -t3, t11, t31; t2, -t1, t10, t30; t43, -t44, -t23, t34; 0, 0, 0, 1; t26 * t11 + t4 * t28, t28 * t11 - t4 * t26, t3, t4 * pkin(4) + t3 * pkin(6) + t31; t26 * t10 + t2 * t28, t28 * t10 - t2 * t26, t1, t2 * pkin(4) + t1 * pkin(6) + t30; -t23 * t26 + t28 * t43, -t23 * t28 - t26 * t43, t44, (pkin(4) * t15 + pkin(6) * t14) * t21 + t34; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
