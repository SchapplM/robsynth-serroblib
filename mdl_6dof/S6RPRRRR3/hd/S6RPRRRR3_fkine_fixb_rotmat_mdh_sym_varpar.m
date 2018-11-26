% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR3
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
% Datum: 2018-11-23 16:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:33:30
% EndTime: 2018-11-23 16:33:30
% DurationCPUTime: 0.16s
% Computational Cost: add. (209->60), mult. (127->70), div. (0->0), fcn. (182->12), ass. (0->41)
t29 = -pkin(9) - pkin(8);
t23 = sin(qJ(4));
t46 = t23 * pkin(4);
t26 = cos(qJ(4));
t10 = t26 * pkin(4) + pkin(3);
t20 = qJ(1) + pkin(11);
t11 = sin(t20);
t45 = t11 * t23;
t27 = cos(qJ(3));
t44 = t11 * t27;
t12 = cos(t20);
t43 = t12 * t27;
t22 = qJ(4) + qJ(5);
t14 = sin(t22);
t42 = t14 * t27;
t15 = cos(t22);
t41 = t15 * t27;
t40 = t23 * t27;
t39 = t26 * t27;
t38 = pkin(6) + 0;
t25 = sin(qJ(1));
t37 = t25 * pkin(1) + 0;
t28 = cos(qJ(1));
t36 = t28 * pkin(1) + 0;
t35 = t11 * pkin(2) + t37;
t13 = qJ(2) + t38;
t34 = t12 * pkin(2) + t11 * pkin(7) + t36;
t24 = sin(qJ(3));
t33 = pkin(3) * t27 + pkin(8) * t24;
t1 = pkin(5) * t15 + t10;
t21 = -pkin(10) + t29;
t32 = t1 * t27 - t21 * t24;
t31 = t10 * t27 - t24 * t29;
t30 = -t12 * pkin(7) + t35;
t16 = qJ(6) + t22;
t9 = cos(t16);
t8 = sin(t16);
t4 = t12 * t24;
t3 = t11 * t24;
t2 = pkin(5) * t14 + t46;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t38; 0, 0, 0, 1; t12, -t11, 0, t36; t11, t12, 0, t37; 0, 0, 1, t13; 0, 0, 0, 1; t43, -t4, t11, t34; t44, -t3, -t12, t30; t24, t27, 0, t13; 0, 0, 0, 1; t12 * t39 + t45, t11 * t26 - t12 * t40, t4, t33 * t12 + t34; t11 * t39 - t12 * t23, -t11 * t40 - t12 * t26, t3, t33 * t11 + t30; t24 * t26, -t24 * t23, -t27, t24 * pkin(3) - t27 * pkin(8) + t13; 0, 0, 0, 1; t11 * t14 + t12 * t41, t11 * t15 - t12 * t42, t4, pkin(4) * t45 + t31 * t12 + t34; t11 * t41 - t12 * t14, -t11 * t42 - t12 * t15, t3 (-pkin(7) - t46) * t12 + t31 * t11 + t35; t24 * t15, -t24 * t14, -t27, t24 * t10 + t27 * t29 + t13; 0, 0, 0, 1; t11 * t8 + t9 * t43, t11 * t9 - t8 * t43, t4, t11 * t2 + t32 * t12 + t34; -t12 * t8 + t9 * t44, -t12 * t9 - t8 * t44, t3 (-pkin(7) - t2) * t12 + t32 * t11 + t35; t24 * t9, -t24 * t8, -t27, t24 * t1 + t27 * t21 + t13; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
