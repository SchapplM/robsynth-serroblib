% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2018-11-23 16:13
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:13:25
% EndTime: 2018-11-23 16:13:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (185->51), mult. (173->47), div. (0->0), fcn. (242->8), ass. (0->38)
t26 = pkin(9) + qJ(3);
t23 = sin(t26);
t31 = sin(qJ(4));
t15 = t23 * t31;
t33 = cos(qJ(4));
t16 = t23 * t33;
t51 = pkin(4) * t16 + qJ(5) * t15;
t32 = sin(qJ(1));
t17 = t32 * t23;
t24 = cos(t26);
t50 = t32 * t24;
t49 = t32 * t31;
t48 = t32 * t33;
t34 = cos(qJ(1));
t18 = t34 * t23;
t47 = t34 * t24;
t46 = t34 * t31;
t45 = t34 * t33;
t44 = qJ(6) * t23;
t27 = pkin(6) + 0;
t28 = sin(pkin(9));
t43 = t28 * pkin(2) + t27;
t29 = cos(pkin(9));
t20 = t29 * pkin(2) + pkin(1);
t30 = -pkin(7) - qJ(2);
t42 = t32 * t20 + t34 * t30 + 0;
t41 = t23 * pkin(3) + t43;
t40 = t34 * t20 - t32 * t30 + 0;
t39 = pkin(3) * t50 + pkin(8) * t17 + t42;
t38 = -t24 * pkin(8) + t41;
t37 = pkin(3) * t47 + pkin(8) * t18 + t40;
t3 = t24 * t49 + t45;
t4 = t24 * t48 - t46;
t36 = t4 * pkin(4) + t3 * qJ(5) + t39;
t5 = t24 * t46 - t48;
t6 = t24 * t45 + t49;
t35 = t6 * pkin(4) + t5 * qJ(5) + t37;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t32, 0, 0; t32, t34, 0, 0; 0, 0, 1, t27; 0, 0, 0, 1; t34 * t29, -t34 * t28, t32, t34 * pkin(1) + t32 * qJ(2) + 0; t32 * t29, -t32 * t28, -t34, t32 * pkin(1) - t34 * qJ(2) + 0; t28, t29, 0, t27; 0, 0, 0, 1; t47, -t18, t32, t40; t50, -t17, -t34, t42; t23, t24, 0, t43; 0, 0, 0, 1; t6, -t5, t18, t37; t4, -t3, t17, t39; t16, -t15, -t24, t38; 0, 0, 0, 1; t6, t18, t5, t35; t4, t17, t3, t36; t16, -t24, t15, t38 + t51; 0, 0, 0, 1; t6, t5, -t18, t6 * pkin(5) - t34 * t44 + t35; t4, t3, -t17, t4 * pkin(5) - t32 * t44 + t36; t16, t15, t24, pkin(5) * t16 + (-pkin(8) + qJ(6)) * t24 + t41 + t51; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
