% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:36
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:35:38
% EndTime: 2018-11-23 16:35:38
% DurationCPUTime: 0.15s
% Computational Cost: add. (199->66), mult. (137->74), div. (0->0), fcn. (196->12), ass. (0->44)
t30 = -pkin(9) - pkin(8);
t26 = sin(qJ(4));
t48 = t26 * pkin(4);
t28 = cos(qJ(4));
t11 = t28 * pkin(4) + pkin(3);
t19 = pkin(11) + qJ(3);
t13 = cos(t19);
t27 = sin(qJ(1));
t47 = t27 * t13;
t22 = qJ(4) + qJ(5);
t14 = sin(t22);
t46 = t27 * t14;
t15 = cos(t22);
t45 = t27 * t15;
t44 = t27 * t26;
t43 = t27 * t28;
t29 = cos(qJ(1));
t42 = t29 * t13;
t41 = t29 * t14;
t40 = t29 * t15;
t39 = t29 * t26;
t38 = t29 * t28;
t20 = pkin(6) + 0;
t24 = cos(pkin(11));
t7 = t24 * pkin(2) + pkin(1);
t37 = t29 * t7 + 0;
t25 = -pkin(7) - qJ(2);
t36 = t29 * t25 + t27 * t7 + 0;
t23 = sin(pkin(11));
t35 = t23 * pkin(2) + t20;
t12 = sin(t19);
t34 = pkin(3) * t13 + pkin(8) * t12;
t1 = pkin(5) * t15 + t11;
t21 = -pkin(10) + t30;
t33 = t1 * t13 - t12 * t21;
t32 = t11 * t13 - t12 * t30;
t31 = -t27 * t25 + t37;
t17 = qJ(6) + t22;
t10 = cos(t17);
t9 = sin(t17);
t6 = t29 * t12;
t5 = t27 * t12;
t2 = pkin(5) * t14 + t48;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t27, 0, 0; t27, t29, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t29 * t24, -t29 * t23, t27, t29 * pkin(1) + t27 * qJ(2) + 0; t27 * t24, -t27 * t23, -t29, t27 * pkin(1) - t29 * qJ(2) + 0; t23, t24, 0, t20; 0, 0, 0, 1; t42, -t6, t27, t31; t47, -t5, -t29, t36; t12, t13, 0, t35; 0, 0, 0, 1; t13 * t38 + t44, -t13 * t39 + t43, t6, t29 * t34 + t31; t13 * t43 - t39, -t13 * t44 - t38, t5, t27 * t34 + t36; t12 * t28, -t12 * t26, -t13, t12 * pkin(3) - t13 * pkin(8) + t35; 0, 0, 0, 1; t13 * t40 + t46, -t13 * t41 + t45, t6, t32 * t29 + (-t25 + t48) * t27 + t37; t13 * t45 - t41, -t13 * t46 - t40, t5, -pkin(4) * t39 + t27 * t32 + t36; t12 * t15, -t12 * t14, -t13, t12 * t11 + t13 * t30 + t35; 0, 0, 0, 1; t10 * t42 + t27 * t9, t27 * t10 - t9 * t42, t6, t33 * t29 + (t2 - t25) * t27 + t37; t10 * t47 - t29 * t9, -t29 * t10 - t9 * t47, t5, -t29 * t2 + t27 * t33 + t36; t12 * t10, -t12 * t9, -t13, t12 * t1 + t13 * t21 + t35; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
