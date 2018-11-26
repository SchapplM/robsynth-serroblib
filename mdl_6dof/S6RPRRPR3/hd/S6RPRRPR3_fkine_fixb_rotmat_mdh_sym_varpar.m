% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2018-11-23 16:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:16:42
% EndTime: 2018-11-23 16:16:42
% DurationCPUTime: 0.14s
% Computational Cost: add. (213->51), mult. (189->54), div. (0->0), fcn. (264->10), ass. (0->37)
t29 = sin(qJ(3));
t32 = cos(qJ(4));
t19 = t29 * t32;
t28 = sin(qJ(4));
t47 = t29 * t28;
t51 = pkin(4) * t19 + qJ(5) * t47;
t26 = qJ(1) + pkin(10);
t20 = sin(t26);
t11 = t20 * t29;
t33 = cos(qJ(3));
t50 = t20 * t33;
t21 = cos(t26);
t13 = t21 * t29;
t49 = t21 * t33;
t48 = t28 * t33;
t46 = t32 * t33;
t45 = pkin(6) + 0;
t30 = sin(qJ(1));
t44 = t30 * pkin(1) + 0;
t34 = cos(qJ(1));
t43 = t34 * pkin(1) + 0;
t22 = qJ(2) + t45;
t42 = t29 * pkin(3) + t22;
t41 = t21 * pkin(2) + t20 * pkin(7) + t43;
t40 = t20 * pkin(2) - t21 * pkin(7) + t44;
t39 = pkin(3) * t49 + pkin(8) * t13 + t41;
t38 = -t33 * pkin(8) + t42;
t37 = pkin(3) * t50 + pkin(8) * t11 + t40;
t5 = -t20 * t32 + t21 * t48;
t6 = t20 * t28 + t21 * t46;
t36 = t6 * pkin(4) + t5 * qJ(5) + t39;
t3 = t20 * t48 + t21 * t32;
t4 = t20 * t46 - t21 * t28;
t35 = t4 * pkin(4) + t3 * qJ(5) + t37;
t31 = cos(qJ(6));
t27 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t30, 0, 0; t30, t34, 0, 0; 0, 0, 1, t45; 0, 0, 0, 1; t21, -t20, 0, t43; t20, t21, 0, t44; 0, 0, 1, t22; 0, 0, 0, 1; t49, -t13, t20, t41; t50, -t11, -t21, t40; t29, t33, 0, t22; 0, 0, 0, 1; t6, -t5, t13, t39; t4, -t3, t11, t37; t19, -t47, -t33, t38; 0, 0, 0, 1; t6, t13, t5, t36; t4, t11, t3, t35; t19, -t33, t47, t38 + t51; 0, 0, 0, 1; t5 * t27 + t6 * t31, -t6 * t27 + t5 * t31, -t13, t6 * pkin(5) - pkin(9) * t13 + t36; t3 * t27 + t4 * t31, -t4 * t27 + t3 * t31, -t11, t4 * pkin(5) - pkin(9) * t11 + t35; (t27 * t28 + t31 * t32) * t29 (-t27 * t32 + t28 * t31) * t29, t33, pkin(5) * t19 + (-pkin(8) + pkin(9)) * t33 + t42 + t51; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
