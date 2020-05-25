% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRP6
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

function [T_c_mdh, Tc_stack] = S6PRPRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:02:24
% EndTime: 2018-11-23 15:02:24
% DurationCPUTime: 0.16s
% Computational Cost: add. (550->70), mult. (613->69), div. (0->0), fcn. (686->14), ass. (0->53)
t41 = sin(pkin(10));
t42 = sin(pkin(6));
t34 = t41 * t42;
t43 = cos(pkin(10));
t73 = t43 * t42;
t72 = pkin(6) - qJ(2);
t71 = pkin(6) + qJ(2);
t70 = pkin(7) * t73;
t69 = t41 * pkin(1) + 0;
t68 = qJ(1) + 0;
t67 = t43 * pkin(1) + pkin(7) * t34 + 0;
t44 = cos(pkin(6));
t66 = t44 * pkin(7) + t68;
t65 = cos(t71);
t64 = sin(t72);
t63 = cos(t72) / 0.2e1;
t62 = sin(t71) / 0.2e1;
t61 = t63 + t65 / 0.2e1;
t47 = sin(qJ(2));
t19 = t41 * t47 - t43 * t61;
t50 = cos(qJ(2));
t57 = t62 - t64 / 0.2e1;
t20 = t41 * t50 + t43 * t57;
t60 = t20 * pkin(2) + t19 * qJ(3) + t69;
t21 = t41 * t61 + t43 * t47;
t22 = -t41 * t57 + t43 * t50;
t59 = t22 * pkin(2) + t21 * qJ(3) + t67;
t30 = t62 + t64 / 0.2e1;
t31 = t63 - t65 / 0.2e1;
t58 = t31 * pkin(2) - t30 * qJ(3) + t66;
t56 = pkin(3) * t34 + t22 * pkin(8) + t59;
t55 = t44 * pkin(3) + t31 * pkin(8) + t58;
t54 = t20 * pkin(8) + (-pkin(3) - pkin(7)) * t73 + t60;
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t10 = t21 * t46 + t49 * t34;
t9 = -t21 * t49 + t46 * t34;
t53 = t10 * pkin(4) + t9 * pkin(9) + t56;
t23 = t30 * t49 + t44 * t46;
t24 = -t30 * t46 + t44 * t49;
t52 = t24 * pkin(4) + t23 * pkin(9) + t55;
t11 = t19 * t49 + t46 * t73;
t12 = t19 * t46 - t49 * t73;
t51 = t12 * pkin(4) - t11 * pkin(9) + t54;
t48 = cos(qJ(5));
t45 = sin(qJ(5));
t6 = t24 * t48 + t31 * t45;
t5 = t24 * t45 - t31 * t48;
t4 = t12 * t48 + t20 * t45;
t3 = t12 * t45 - t20 * t48;
t2 = t10 * t48 + t22 * t45;
t1 = t10 * t45 - t22 * t48;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t43, -t41, 0, 0; t41, t43, 0, 0; 0, 0, 1, t68; 0, 0, 0, 1; t22, -t21, t34, t67; t20, -t19, -t73, t69 - t70; t31, t30, t44, t66; 0, 0, 0, 1; t34, -t22, t21, t59; -t73, -t20, t19, t60 - t70; t44, -t31, -t30, t58; 0, 0, 0, 1; t10, -t9, t22, t56; t12, t11, t20, t54; t24, -t23, t31, t55; 0, 0, 0, 1; t2, -t1, t9, t53; t4, -t3, -t11, t51; t6, -t5, t23, t52; 0, 0, 0, 1; t2, t9, t1, t2 * pkin(5) + t1 * qJ(6) + t53; t4, -t11, t3, t4 * pkin(5) + t3 * qJ(6) + t51; t6, t23, t5, t6 * pkin(5) + t5 * qJ(6) + t52; 0, 0, 0, 1;];
T_ges = t7;
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
