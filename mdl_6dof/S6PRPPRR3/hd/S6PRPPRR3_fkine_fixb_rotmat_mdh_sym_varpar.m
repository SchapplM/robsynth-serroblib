% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 14:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRPPRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:54:27
% EndTime: 2018-11-23 14:54:27
% DurationCPUTime: 0.17s
% Computational Cost: add. (688->72), mult. (741->82), div. (0->0), fcn. (818->16), ass. (0->55)
t70 = pkin(6) + qJ(2);
t75 = sin(t70) / 0.2e1;
t43 = sin(pkin(10));
t44 = sin(pkin(6));
t36 = t43 * t44;
t45 = cos(pkin(10));
t74 = t45 * t44;
t73 = qJ(4) * t44;
t72 = cos(pkin(11));
t71 = pkin(6) - qJ(2);
t69 = qJ(1) + 0;
t68 = t45 * pkin(1) + pkin(7) * t36 + 0;
t46 = cos(pkin(6));
t67 = t46 * pkin(7) + t69;
t66 = cos(t70);
t65 = sin(t71);
t64 = cos(t71) / 0.2e1;
t63 = t43 * pkin(1) - pkin(7) * t74 + 0;
t49 = sin(qJ(2));
t59 = t64 + t66 / 0.2e1;
t24 = t43 * t59 + t45 * t49;
t31 = t75 - t65 / 0.2e1;
t52 = cos(qJ(2));
t25 = -t43 * t31 + t45 * t52;
t62 = t25 * pkin(2) + t24 * qJ(3) + t68;
t30 = t75 + t65 / 0.2e1;
t32 = t64 - t66 / 0.2e1;
t61 = t32 * pkin(2) - t30 * qJ(3) + t67;
t22 = t43 * t49 - t45 * t59;
t23 = t45 * t31 + t43 * t52;
t60 = t23 * pkin(2) + t22 * qJ(3) + t63;
t58 = t32 * pkin(3) - t46 * qJ(4) + t61;
t57 = t23 * pkin(3) + t45 * t73 + t60;
t56 = t25 * pkin(3) - t43 * t73 + t62;
t42 = sin(pkin(11));
t14 = t30 * t72 + t32 * t42;
t15 = -t30 * t42 + t32 * t72;
t55 = t15 * pkin(4) + t14 * pkin(8) + t58;
t7 = -t22 * t72 + t23 * t42;
t8 = t22 * t42 + t23 * t72;
t54 = t8 * pkin(4) + t7 * pkin(8) + t57;
t10 = t24 * t42 + t25 * t72;
t9 = -t24 * t72 + t25 * t42;
t53 = t10 * pkin(4) + t9 * pkin(8) + t56;
t51 = cos(qJ(5));
t50 = cos(qJ(6));
t48 = sin(qJ(5));
t47 = sin(qJ(6));
t12 = t15 * t51 - t46 * t48;
t11 = t15 * t48 + t46 * t51;
t4 = t10 * t51 - t48 * t36;
t3 = t10 * t48 + t51 * t36;
t2 = t48 * t74 + t8 * t51;
t1 = t8 * t48 - t51 * t74;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t45, -t43, 0, 0; t43, t45, 0, 0; 0, 0, 1, t69; 0, 0, 0, 1; t25, -t24, t36, t68; t23, -t22, -t74, t63; t32, t30, t46, t67; 0, 0, 0, 1; t25, t36, t24, t62; t23, -t74, t22, t60; t32, t46, -t30, t61; 0, 0, 0, 1; t10, -t9, -t36, t56; t8, -t7, t74, t57; t15, -t14, -t46, t58; 0, 0, 0, 1; t4, -t3, t9, t53; t2, -t1, t7, t54; t12, -t11, t14, t55; 0, 0, 0, 1; t4 * t50 + t9 * t47, -t4 * t47 + t9 * t50, t3, t4 * pkin(5) + t3 * pkin(9) + t53; t2 * t50 + t7 * t47, -t2 * t47 + t7 * t50, t1, t2 * pkin(5) + t1 * pkin(9) + t54; t12 * t50 + t14 * t47, -t12 * t47 + t14 * t50, t11, t12 * pkin(5) + t11 * pkin(9) + t55; 0, 0, 0, 1;];
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
