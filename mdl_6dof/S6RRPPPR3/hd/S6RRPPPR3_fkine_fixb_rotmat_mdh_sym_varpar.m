% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2018-11-23 16:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPPPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:43:24
% EndTime: 2018-11-23 16:43:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (116->59), mult. (152->56), div. (0->0), fcn. (206->8), ass. (0->35)
t26 = sin(qJ(1));
t25 = sin(qJ(2));
t42 = qJ(3) * t25;
t27 = cos(qJ(2));
t9 = t26 * t27;
t45 = pkin(2) * t9 + t26 * t42;
t28 = cos(qJ(1));
t44 = pkin(3) * t9 + t28 * qJ(4);
t22 = sin(pkin(9));
t43 = pkin(5) * t22;
t8 = t26 * t25;
t10 = t28 * t25;
t11 = t28 * t27;
t21 = pkin(6) + 0;
t41 = t26 * pkin(1) + 0;
t40 = t25 * pkin(2) + t21;
t39 = t28 * pkin(1) + t26 * pkin(7) + 0;
t15 = t25 * pkin(3);
t38 = t15 + t40;
t24 = -pkin(8) - qJ(5);
t23 = cos(pkin(9));
t7 = t23 * pkin(5) + pkin(4);
t37 = -t24 * t27 + t25 * t7;
t36 = pkin(4) * t25 + qJ(5) * t27;
t35 = -t28 * pkin(7) + t41;
t34 = pkin(2) * t11 + t28 * t42 + t39;
t33 = pkin(3) * t11 + t34;
t32 = -t27 * qJ(3) + t40;
t31 = t35 + t45;
t30 = t31 + t44;
t29 = -t26 * qJ(4) + t33;
t20 = pkin(9) + qJ(6);
t13 = cos(t20);
t12 = sin(t20);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t26, 0, 0; t26, t28, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t11, -t10, t26, t39; t9, -t8, -t28, t35; t25, t27, 0, t21; 0, 0, 0, 1; t11, t26, t10, t34; t9, -t28, t8, t31; t25, 0, -t27, t32; 0, 0, 0, 1; t10, -t11, -t26, t29; t8, -t9, t28, t30; -t27, -t25, 0, t15 + t32; 0, 0, 0, 1; t23 * t10 - t26 * t22, -t22 * t10 - t26 * t23, t11, t36 * t28 + t29; t28 * t22 + t23 * t8, -t22 * t8 + t28 * t23, t9, t36 * t26 + t30; -t27 * t23, t27 * t22, t25, t25 * qJ(5) + (-pkin(4) - qJ(3)) * t27 + t38; 0, 0, 0, 1; t13 * t10 - t26 * t12, -t12 * t10 - t26 * t13, t11, t37 * t28 + (-qJ(4) - t43) * t26 + t33; t28 * t12 + t13 * t8, -t12 * t8 + t28 * t13, t9 (-pkin(7) + t43) * t28 + t37 * t26 + t41 + t44 + t45; -t27 * t13, t27 * t12, t25, -t25 * t24 + (-qJ(3) - t7) * t27 + t38; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
