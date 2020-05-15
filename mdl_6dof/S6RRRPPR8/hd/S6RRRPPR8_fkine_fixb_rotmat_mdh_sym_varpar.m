% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2018-11-23 17:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:38:23
% EndTime: 2018-11-23 17:38:24
% DurationCPUTime: 0.18s
% Computational Cost: add. (552->74), mult. (621->70), div. (0->0), fcn. (697->14), ass. (0->55)
t70 = cos(qJ(3));
t33 = sin(pkin(6));
t37 = sin(qJ(1));
t69 = t37 * t33;
t40 = cos(qJ(1));
t68 = t40 * t33;
t60 = pkin(6) + qJ(2);
t50 = sin(t60) / 0.2e1;
t61 = pkin(6) - qJ(2);
t54 = sin(t61);
t22 = t50 - t54 / 0.2e1;
t39 = cos(qJ(2));
t16 = t40 * t22 + t37 * t39;
t35 = sin(qJ(3));
t58 = t33 * t70;
t5 = t16 * t35 + t40 * t58;
t67 = t5 * qJ(4);
t18 = -t37 * t22 + t40 * t39;
t7 = t18 * t35 - t37 * t58;
t66 = t7 * qJ(4);
t65 = pkin(5) + qJ(4);
t64 = pkin(9) - qJ(5);
t51 = cos(t61) / 0.2e1;
t55 = cos(t60);
t23 = t51 - t55 / 0.2e1;
t62 = cos(pkin(6));
t13 = t23 * t35 - t62 * t70;
t63 = t13 * qJ(4);
t59 = pkin(7) + 0;
t57 = t62 * pkin(8) + t59;
t56 = t40 * pkin(1) + pkin(8) * t69 + 0;
t53 = t23 * pkin(2) + t57;
t52 = t18 * pkin(2) + t56;
t49 = t37 * pkin(1) - pkin(8) * t68 + 0;
t21 = t50 + t54 / 0.2e1;
t48 = -t21 * pkin(9) + t53;
t47 = t16 * pkin(2) + t49;
t36 = sin(qJ(2));
t44 = t51 + t55 / 0.2e1;
t17 = t40 * t36 + t37 * t44;
t46 = t17 * pkin(9) + t52;
t15 = t37 * t36 - t40 * t44;
t45 = t15 * pkin(9) + t47;
t14 = t23 * t70 + t62 * t35;
t10 = t14 * pkin(3);
t43 = t14 * pkin(4) - t64 * t21 + t10 + t53;
t8 = t18 * t70 + t35 * t69;
t4 = t8 * pkin(3);
t42 = t8 * pkin(4) + t64 * t17 + t4 + t52;
t6 = t16 * t70 - t35 * t68;
t2 = t6 * pkin(3);
t41 = t6 * pkin(4) + t64 * t15 + t2 + t47;
t38 = cos(qJ(6));
t34 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t40, -t37, 0, 0; t37, t40, 0, 0; 0, 0, 1, t59; 0, 0, 0, 1; t18, -t17, t69, t56; t16, -t15, -t68, t49; t23, t21, t62, t57; 0, 0, 0, 1; t8, -t7, t17, t46; t6, -t5, t15, t45; t14, -t13, -t21, t48; 0, 0, 0, 1; t8, t17, t7, t4 + t46 + t66; t6, t15, t5, t2 + t45 + t67; t14, -t21, t13, t10 + t48 + t63; 0, 0, 0, 1; t7, -t8, -t17, t42 + t66; t5, -t6, -t15, t41 + t67; t13, -t14, t21, t43 + t63; 0, 0, 0, 1; -t17 * t34 + t7 * t38, -t17 * t38 - t7 * t34, t8, t8 * pkin(10) + t65 * t7 + t42; -t15 * t34 + t5 * t38, -t15 * t38 - t5 * t34, t6, t6 * pkin(10) + t65 * t5 + t41; t13 * t38 + t21 * t34, -t13 * t34 + t21 * t38, t14, t14 * pkin(10) + t65 * t13 + t43; 0, 0, 0, 1;];
T_ges = t1;
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
