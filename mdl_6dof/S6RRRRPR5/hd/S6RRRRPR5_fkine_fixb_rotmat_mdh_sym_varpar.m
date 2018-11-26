% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:14:35
% EndTime: 2018-11-23 18:14:35
% DurationCPUTime: 0.13s
% Computational Cost: add. (201->57), mult. (199->60), div. (0->0), fcn. (278->10), ass. (0->39)
t26 = qJ(2) + qJ(3);
t22 = sin(t26);
t32 = cos(qJ(4));
t15 = t22 * t32;
t28 = sin(qJ(4));
t51 = t22 * t28;
t52 = pkin(4) * t15 + qJ(5) * t51;
t30 = sin(qJ(1));
t16 = t30 * t22;
t23 = cos(t26);
t50 = t30 * t23;
t49 = t30 * t28;
t48 = t30 * t32;
t34 = cos(qJ(1));
t17 = t34 * t22;
t47 = t34 * t23;
t46 = t34 * t28;
t45 = t34 * t32;
t25 = pkin(6) + 0;
t29 = sin(qJ(2));
t44 = t29 * pkin(2) + t25;
t33 = cos(qJ(2));
t20 = t33 * pkin(2) + pkin(1);
t35 = -pkin(8) - pkin(7);
t43 = t30 * t20 + t34 * t35 + 0;
t42 = t22 * pkin(3) + t44;
t41 = t34 * t20 - t30 * t35 + 0;
t40 = pkin(3) * t50 + pkin(9) * t16 + t43;
t39 = -t23 * pkin(9) + t42;
t38 = pkin(3) * t47 + pkin(9) * t17 + t41;
t3 = t23 * t49 + t45;
t4 = t23 * t48 - t46;
t37 = t4 * pkin(4) + t3 * qJ(5) + t40;
t5 = t23 * t46 - t48;
t6 = t23 * t45 + t49;
t36 = t6 * pkin(4) + t5 * qJ(5) + t38;
t31 = cos(qJ(6));
t27 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t30, 0, 0; t30, t34, 0, 0; 0, 0, 1, t25; 0, 0, 0, 1; t34 * t33, -t34 * t29, t30, t34 * pkin(1) + t30 * pkin(7) + 0; t30 * t33, -t30 * t29, -t34, t30 * pkin(1) - t34 * pkin(7) + 0; t29, t33, 0, t25; 0, 0, 0, 1; t47, -t17, t30, t41; t50, -t16, -t34, t43; t22, t23, 0, t44; 0, 0, 0, 1; t6, -t5, t17, t38; t4, -t3, t16, t40; t15, -t51, -t23, t39; 0, 0, 0, 1; t6, t17, t5, t36; t4, t16, t3, t37; t15, -t23, t51, t39 + t52; 0, 0, 0, 1; t5 * t27 + t6 * t31, -t6 * t27 + t5 * t31, -t17, t6 * pkin(5) - pkin(10) * t17 + t36; t3 * t27 + t4 * t31, -t4 * t27 + t3 * t31, -t16, t4 * pkin(5) - pkin(10) * t16 + t37; (t27 * t28 + t31 * t32) * t22 (-t27 * t32 + t28 * t31) * t22, t23, pkin(5) * t15 + (-pkin(9) + pkin(10)) * t23 + t42 + t52; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
