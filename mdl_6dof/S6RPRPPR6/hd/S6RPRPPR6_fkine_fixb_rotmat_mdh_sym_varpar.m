% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2018-11-23 15:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRPPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:54:36
% EndTime: 2018-11-23 15:54:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (141->56), mult. (110->56), div. (0->0), fcn. (160->10), ass. (0->40)
t20 = sin(qJ(1));
t13 = qJ(3) + pkin(9);
t5 = sin(t13);
t44 = t20 * t5;
t12 = pkin(10) + qJ(6);
t6 = cos(t12);
t43 = t20 * t6;
t7 = cos(t13);
t42 = t20 * t7;
t22 = cos(qJ(1));
t41 = t22 * t5;
t40 = t22 * t6;
t15 = sin(pkin(10));
t39 = t20 * t15;
t16 = cos(pkin(10));
t38 = t20 * t16;
t19 = sin(qJ(3));
t37 = t20 * t19;
t36 = t22 * t15;
t35 = t22 * t16;
t14 = pkin(6) + 0;
t34 = t20 * pkin(1) + 0;
t33 = pkin(2) + t14;
t17 = -qJ(4) - pkin(7);
t32 = pkin(5) * t15 - t17;
t31 = t22 * pkin(1) + t20 * qJ(2) + 0;
t30 = -pkin(3) * t19 - qJ(2);
t21 = cos(qJ(3));
t29 = t21 * pkin(3) + t33;
t28 = pkin(3) * t37 + t31;
t18 = -pkin(8) - qJ(5);
t3 = t16 * pkin(5) + pkin(4);
t27 = t18 * t7 + t3 * t5;
t26 = pkin(4) * t5 - qJ(5) * t7;
t25 = -t20 * t17 + t34;
t24 = -t22 * qJ(2) + t34;
t23 = -t22 * t17 + t28;
t4 = sin(t12);
t1 = t22 * t7;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t20, 0, 0; t20, t22, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; 0, -t22, t20, t31; 0, -t20, -t22, t24; 1, 0, 0, t14; 0, 0, 0, 1; t37, t20 * t21, t22, t22 * pkin(7) + t31; -t22 * t19, -t22 * t21, t20, t20 * pkin(7) + t24; t21, -t19, 0, t33; 0, 0, 0, 1; t44, t42, t22, t23; -t41, -t1, t20, t30 * t22 + t25; t7, -t5, 0, t29; 0, 0, 0, 1; t5 * t38 + t36, -t5 * t39 + t35, -t42, t26 * t20 + t23; -t5 * t35 + t39, t5 * t36 + t38, t1 (-t26 + t30) * t22 + t25; t7 * t16, -t7 * t15, t5, t7 * pkin(4) + t5 * qJ(5) + t29; 0, 0, 0, 1; t22 * t4 + t5 * t43, -t4 * t44 + t40, -t42, t27 * t20 + t32 * t22 + t28; t20 * t4 - t5 * t40, t4 * t41 + t43, t1, t32 * t20 + (-t27 + t30) * t22 + t34; t7 * t6, -t7 * t4, t5, -t5 * t18 + t7 * t3 + t29; 0, 0, 0, 1;];
T_ges = t2;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
