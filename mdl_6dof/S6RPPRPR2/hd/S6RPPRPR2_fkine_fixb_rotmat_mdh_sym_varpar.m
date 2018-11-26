% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2018-11-23 15:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPPRPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:39:58
% EndTime: 2018-11-23 15:39:58
% DurationCPUTime: 0.14s
% Computational Cost: add. (184->49), mult. (90->42), div. (0->0), fcn. (135->10), ass. (0->35)
t21 = qJ(1) + pkin(9);
t15 = cos(t21);
t20 = pkin(10) + qJ(4);
t12 = sin(t20);
t38 = qJ(5) * t12;
t14 = cos(t20);
t8 = t15 * t14;
t45 = pkin(4) * t8 + t15 * t38;
t13 = sin(t21);
t44 = t13 * t12;
t7 = t13 * t14;
t25 = sin(qJ(6));
t43 = t13 * t25;
t27 = cos(qJ(6));
t42 = t13 * t27;
t41 = t15 * t12;
t40 = t15 * t25;
t39 = t15 * t27;
t37 = pkin(6) + 0;
t26 = sin(qJ(1));
t36 = t26 * pkin(1) + 0;
t28 = cos(qJ(1));
t35 = t28 * pkin(1) + 0;
t23 = cos(pkin(10));
t11 = t23 * pkin(3) + pkin(2);
t34 = t15 * t11 + t35;
t16 = qJ(2) + t37;
t24 = -pkin(7) - qJ(3);
t33 = t13 * t11 + t15 * t24 + t36;
t22 = sin(pkin(10));
t32 = t22 * pkin(3) + t16;
t31 = pkin(4) * t7 + t13 * t38 + t33;
t30 = -t13 * t24 + t34;
t29 = t12 * pkin(4) - t14 * qJ(5) + t32;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t26, 0, 0; t26, t28, 0, 0; 0, 0, 1, t37; 0, 0, 0, 1; t15, -t13, 0, t35; t13, t15, 0, t36; 0, 0, 1, t16; 0, 0, 0, 1; t15 * t23, -t15 * t22, t13, t15 * pkin(2) + t13 * qJ(3) + t35; t13 * t23, -t13 * t22, -t15, t13 * pkin(2) - t15 * qJ(3) + t36; t22, t23, 0, t16; 0, 0, 0, 1; t8, -t41, t13, t30; t7, -t44, -t15, t33; t12, t14, 0, t32; 0, 0, 0, 1; t13, -t8, t41, t30 + t45; -t15, -t7, t44, t31; 0, -t12, -t14, t29; 0, 0, 0, 1; t12 * t40 + t42, t12 * t39 - t43, t8, pkin(8) * t8 + (pkin(5) - t24) * t13 + t34 + t45; t12 * t43 - t39, t12 * t42 + t40, t7, -t15 * pkin(5) + pkin(8) * t7 + t31; -t14 * t25, -t14 * t27, t12, t12 * pkin(8) + t29; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
