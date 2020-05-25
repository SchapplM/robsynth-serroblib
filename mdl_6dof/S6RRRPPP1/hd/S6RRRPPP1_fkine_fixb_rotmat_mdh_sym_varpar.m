% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2018-11-23 17:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPPP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:31:49
% EndTime: 2018-11-23 17:31:49
% DurationCPUTime: 0.18s
% Computational Cost: add. (445->70), mult. (658->80), div. (0->0), fcn. (786->14), ass. (0->57)
t44 = sin(qJ(3));
t45 = sin(qJ(2));
t79 = t45 * t44;
t47 = cos(qJ(3));
t78 = t45 * t47;
t46 = sin(qJ(1));
t77 = t46 * t45;
t48 = cos(qJ(2));
t76 = t46 * t48;
t49 = cos(qJ(1));
t75 = t49 * t45;
t74 = t49 * t48;
t41 = sin(pkin(6));
t73 = qJ(4) * t41;
t43 = cos(pkin(6));
t72 = qJ(4) * t43;
t39 = pkin(7) + 0;
t71 = t41 * t79;
t70 = pkin(6) - pkin(10);
t69 = pkin(6) + pkin(10);
t68 = t45 * t72;
t67 = t45 * pkin(2) + t39;
t66 = t49 * pkin(1) + t46 * pkin(8) + 0;
t65 = cos(t69);
t64 = sin(t70);
t63 = t46 * pkin(1) - t49 * pkin(8) + 0;
t62 = cos(t70) / 0.2e1;
t61 = sin(t69) / 0.2e1;
t60 = pkin(2) * t74 + pkin(9) * t75 + t66;
t59 = pkin(2) * t76 + pkin(9) * t77 + t63;
t58 = t62 + t65 / 0.2e1;
t57 = t61 + t64 / 0.2e1;
t56 = qJ(4) * t71 + pkin(3) * t78 + (-pkin(9) - t72) * t48 + t67;
t23 = -t44 * t74 + t46 * t47;
t24 = t46 * t44 + t47 * t74;
t55 = t24 * pkin(3) - t23 * t73 + t49 * t68 + t60;
t54 = t45 * t57;
t21 = -t44 * t76 - t49 * t47;
t22 = -t49 * t44 + t47 * t76;
t53 = t22 * pkin(3) - t21 * t73 + t46 * t68 + t59;
t40 = sin(pkin(10));
t8 = t48 * t57 + (t47 * t40 + t44 * t58) * t45;
t25 = t61 - t64 / 0.2e1;
t26 = t62 - t65 / 0.2e1;
t42 = cos(pkin(10));
t9 = -t25 * t79 - t48 * t26 + t42 * t78;
t52 = t9 * pkin(4) + t8 * qJ(5) + t56;
t5 = -t23 * t58 + t24 * t40 - t49 * t54;
t6 = t23 * t25 + t24 * t42 + t26 * t75;
t51 = t6 * pkin(4) + t5 * qJ(5) + t55;
t3 = -t21 * t58 + t22 * t40 - t46 * t54;
t4 = t21 * t25 + t22 * t42 + t26 * t77;
t50 = t4 * pkin(4) + t3 * qJ(5) + t53;
t19 = -t48 * t43 + t71;
t11 = -t23 * t41 + t43 * t75;
t10 = -t21 * t41 + t43 * t77;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t49, -t46, 0, 0; t46, t49, 0, 0; 0, 0, 1, t39; 0, 0, 0, 1; t74, -t75, t46, t66; t76, -t77, -t49, t63; t45, t48, 0, t39; 0, 0, 0, 1; t24, t23, t75, t60; t22, t21, t77, t59; t78, -t79, -t48, -t48 * pkin(9) + t67; 0, 0, 0, 1; t6, -t5, t11, t55; t4, -t3, t10, t53; t9, -t8, t19, t56; 0, 0, 0, 1; t11, -t6, t5, t51; t10, -t4, t3, t50; t19, -t9, t8, t52; 0, 0, 0, 1; t11, t5, t6, t11 * pkin(5) + t6 * qJ(6) + t51; t10, t3, t4, t10 * pkin(5) + t4 * qJ(6) + t50; t19, t8, t9, t19 * pkin(5) + t9 * qJ(6) + t52; 0, 0, 0, 1;];
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
