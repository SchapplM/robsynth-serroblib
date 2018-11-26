% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2018-11-23 15:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRPPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:54:12
% EndTime: 2018-11-23 15:54:12
% DurationCPUTime: 0.15s
% Computational Cost: add. (173->62), mult. (130->61), div. (0->0), fcn. (184->10), ass. (0->41)
t28 = cos(qJ(1));
t19 = pkin(9) + qJ(3);
t14 = sin(t19);
t39 = qJ(4) * t14;
t16 = cos(t19);
t8 = t28 * t16;
t49 = pkin(3) * t8 + t28 * t39;
t21 = sin(pkin(10));
t48 = pkin(5) * t21;
t27 = sin(qJ(1));
t47 = t27 * t14;
t18 = pkin(10) + qJ(6);
t15 = cos(t18);
t46 = t27 * t15;
t7 = t27 * t16;
t45 = t27 * t21;
t23 = cos(pkin(10));
t44 = t27 * t23;
t43 = t28 * t14;
t42 = t28 * t15;
t41 = t28 * t21;
t40 = t28 * t23;
t38 = qJ(5) * t16;
t20 = pkin(6) + 0;
t24 = cos(pkin(9));
t11 = t24 * pkin(2) + pkin(1);
t37 = t28 * t11 + 0;
t22 = sin(pkin(9));
t36 = t22 * pkin(2) + t20;
t26 = -pkin(7) - qJ(2);
t35 = t27 * t11 + t28 * t26 + 0;
t34 = t37 + t49;
t33 = t14 * pkin(3) + t36;
t32 = pkin(3) * t7 + t27 * t39 + t35;
t31 = -t27 * t26 + t37;
t25 = -pkin(8) - qJ(5);
t30 = t14 * t48 - t16 * t25;
t29 = -t16 * qJ(4) + t33;
t13 = sin(t18);
t10 = t23 * pkin(5) + pkin(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t27, 0, 0; t27, t28, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t28 * t24, -t28 * t22, t27, t28 * pkin(1) + t27 * qJ(2) + 0; t27 * t24, -t27 * t22, -t28, t27 * pkin(1) - t28 * qJ(2) + 0; t22, t24, 0, t20; 0, 0, 0, 1; t8, -t43, t27, t31; t7, -t47, -t28, t35; t14, t16, 0, t36; 0, 0, 0, 1; t27, -t8, t43, t31 + t49; -t28, -t7, t47, t32; 0, -t14, -t16, t29; 0, 0, 0, 1; t14 * t41 + t44, t14 * t40 - t45, t8, t28 * t38 + (pkin(4) - t26) * t27 + t34; t14 * t45 - t40, t14 * t44 + t41, t7, -t28 * pkin(4) + t27 * t38 + t32; -t16 * t21, -t16 * t23, t14, t14 * qJ(5) + t29; 0, 0, 0, 1; t13 * t43 + t46, -t27 * t13 + t14 * t42, t8, t30 * t28 + (t10 - t26) * t27 + t34; t13 * t47 - t42, t28 * t13 + t14 * t46, t7, -t28 * t10 + t27 * t30 + t32; -t16 * t13, -t16 * t15, t14, -t14 * t25 + (-qJ(4) - t48) * t16 + t33; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
