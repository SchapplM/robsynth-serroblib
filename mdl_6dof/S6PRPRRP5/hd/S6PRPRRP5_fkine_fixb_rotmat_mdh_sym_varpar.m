% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRPRRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:01:56
% EndTime: 2018-11-23 15:01:56
% DurationCPUTime: 0.18s
% Computational Cost: add. (508->73), mult. (558->74), div. (0->0), fcn. (622->14), ass. (0->56)
t36 = sin(pkin(10));
t37 = sin(pkin(6));
t28 = t36 * t37;
t38 = cos(pkin(10));
t70 = t38 * t37;
t69 = pkin(6) - qJ(2);
t68 = pkin(6) + qJ(2);
t67 = pkin(7) * t70;
t66 = t36 * pkin(1) + 0;
t65 = qJ(1) + 0;
t41 = sin(qJ(5));
t64 = pkin(5) * t41 + pkin(8);
t63 = t38 * pkin(1) + pkin(7) * t28 + 0;
t39 = cos(pkin(6));
t62 = t39 * pkin(7) + t65;
t61 = cos(t68);
t60 = sin(t69);
t59 = cos(t69) / 0.2e1;
t58 = sin(t68) / 0.2e1;
t57 = t59 + t61 / 0.2e1;
t43 = sin(qJ(2));
t14 = t36 * t43 - t38 * t57;
t46 = cos(qJ(2));
t51 = t58 - t60 / 0.2e1;
t15 = t36 * t46 + t38 * t51;
t56 = t15 * pkin(2) + t14 * qJ(3) + t66;
t16 = t36 * t57 + t38 * t43;
t17 = -t36 * t51 + t38 * t46;
t55 = t17 * pkin(2) + t16 * qJ(3) + t63;
t24 = t58 + t60 / 0.2e1;
t25 = t59 - t61 / 0.2e1;
t54 = t25 * pkin(2) - t24 * qJ(3) + t62;
t53 = pkin(3) * t28 + t55;
t52 = t39 * pkin(3) + t54;
t50 = t17 * pkin(8) + t53;
t49 = t25 * pkin(8) + t52;
t48 = (-pkin(3) - pkin(7)) * t70 + t56;
t47 = t15 * pkin(8) + t48;
t45 = cos(qJ(4));
t44 = cos(qJ(5));
t42 = sin(qJ(4));
t40 = -qJ(6) - pkin(9);
t31 = t44 * pkin(5) + pkin(4);
t19 = -t24 * t42 + t39 * t45;
t18 = t24 * t45 + t39 * t42;
t10 = t14 * t42 - t45 * t70;
t9 = t14 * t45 + t42 * t70;
t8 = t16 * t42 + t45 * t28;
t7 = -t16 * t45 + t42 * t28;
t6 = t19 * t44 + t25 * t41;
t5 = -t19 * t41 + t25 * t44;
t4 = t10 * t44 + t15 * t41;
t3 = -t10 * t41 + t15 * t44;
t2 = t17 * t41 + t8 * t44;
t1 = t17 * t44 - t8 * t41;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t38, -t36, 0, 0; t36, t38, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t17, -t16, t28, t63; t15, -t14, -t70, t66 - t67; t25, t24, t39, t62; 0, 0, 0, 1; t28, -t17, t16, t55; -t70, -t15, t14, t56 - t67; t39, -t25, -t24, t54; 0, 0, 0, 1; t8, -t7, t17, t50; t10, t9, t15, t47; t19, -t18, t25, t49; 0, 0, 0, 1; t2, t1, t7, t8 * pkin(4) + t7 * pkin(9) + t50; t4, t3, -t9, t10 * pkin(4) - t9 * pkin(9) + t47; t6, t5, t18, t19 * pkin(4) + t18 * pkin(9) + t49; 0, 0, 0, 1; t2, t1, t7, t64 * t17 + t8 * t31 - t7 * t40 + t53; t4, t3, -t9, t10 * t31 + t64 * t15 + t9 * t40 + t48; t6, t5, t18, -t18 * t40 + t19 * t31 + t64 * t25 + t52; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
