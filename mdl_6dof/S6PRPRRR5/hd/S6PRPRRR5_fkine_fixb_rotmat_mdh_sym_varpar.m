% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:05:52
% EndTime: 2018-11-23 15:05:52
% DurationCPUTime: 0.19s
% Computational Cost: add. (505->79), mult. (495->85), div. (0->0), fcn. (547->16), ass. (0->52)
t40 = sin(pkin(11));
t42 = cos(pkin(11));
t46 = sin(qJ(2));
t69 = pkin(6) - qJ(2);
t60 = cos(t69) / 0.2e1;
t68 = pkin(6) + qJ(2);
t62 = cos(t68);
t53 = t60 + t62 / 0.2e1;
t14 = t40 * t46 - t42 * t53;
t45 = sin(qJ(4));
t73 = t14 * t45;
t16 = t40 * t53 + t42 * t46;
t72 = t16 * t45;
t59 = sin(t68) / 0.2e1;
t61 = sin(t69);
t25 = t59 + t61 / 0.2e1;
t71 = t25 * t45;
t41 = sin(pkin(6));
t29 = t40 * t41;
t70 = t42 * t41;
t67 = pkin(7) * t70;
t66 = t40 * pkin(1) + 0;
t65 = qJ(1) + 0;
t64 = t42 * pkin(1) + pkin(7) * t29 + 0;
t43 = cos(pkin(6));
t63 = t43 * pkin(7) + t65;
t58 = t59 - t61 / 0.2e1;
t49 = cos(qJ(2));
t15 = t40 * t49 + t42 * t58;
t57 = t15 * pkin(2) + t14 * qJ(3) + t66;
t17 = -t40 * t58 + t42 * t49;
t56 = t17 * pkin(2) + t16 * qJ(3) + t64;
t26 = t60 - t62 / 0.2e1;
t55 = t26 * pkin(2) - t25 * qJ(3) + t63;
t48 = cos(qJ(4));
t33 = t48 * pkin(4) + pkin(3);
t50 = -pkin(9) - pkin(8);
t54 = pkin(4) * t72 - t17 * t50 + t33 * t29 + t56;
t52 = -pkin(4) * t71 - t26 * t50 + t43 * t33 + t55;
t51 = -t15 * t50 + pkin(4) * t73 + (-pkin(7) - t33) * t70 + t57;
t47 = cos(qJ(6));
t44 = sin(qJ(6));
t39 = qJ(4) + qJ(5);
t35 = cos(t39);
t34 = sin(t39);
t9 = -t25 * t34 + t43 * t35;
t8 = t25 * t35 + t43 * t34;
t4 = t14 * t34 - t35 * t70;
t3 = t14 * t35 + t34 * t70;
t2 = t16 * t34 + t35 * t29;
t1 = -t16 * t35 + t34 * t29;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t42, -t40, 0, 0; t40, t42, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t17, -t16, t29, t64; t15, -t14, -t70, t66 - t67; t26, t25, t43, t63; 0, 0, 0, 1; t29, -t17, t16, t56; -t70, -t15, t14, t57 - t67; t43, -t26, -t25, t55; 0, 0, 0, 1; t48 * t29 + t72, t16 * t48 - t45 * t29, t17, pkin(3) * t29 + t17 * pkin(8) + t56; -t48 * t70 + t73, t14 * t48 + t45 * t70, t15, t15 * pkin(8) + (-pkin(3) - pkin(7)) * t70 + t57; t43 * t48 - t71, -t25 * t48 - t43 * t45, t26, t43 * pkin(3) + t26 * pkin(8) + t55; 0, 0, 0, 1; t2, -t1, t17, t54; t4, t3, t15, t51; t9, -t8, t26, t52; 0, 0, 0, 1; t17 * t44 + t2 * t47, t17 * t47 - t2 * t44, t1, t2 * pkin(5) + t1 * pkin(10) + t54; t15 * t44 + t4 * t47, t15 * t47 - t4 * t44, -t3, t4 * pkin(5) - t3 * pkin(10) + t51; t26 * t44 + t9 * t47, t26 * t47 - t9 * t44, t8, t9 * pkin(5) + t8 * pkin(10) + t52; 0, 0, 0, 1;];
T_ges = t5;
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
