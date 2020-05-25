% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2018-11-23 15:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRRPRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:13:44
% EndTime: 2018-11-23 15:13:44
% DurationCPUTime: 0.16s
% Computational Cost: add. (611->67), mult. (693->70), div. (0->0), fcn. (780->14), ass. (0->53)
t74 = cos(qJ(3));
t42 = sin(pkin(10));
t43 = sin(pkin(6));
t73 = t42 * t43;
t44 = cos(pkin(10));
t72 = t44 * t43;
t71 = cos(pkin(6));
t70 = pkin(6) - qJ(2);
t69 = pkin(6) + qJ(2);
t68 = qJ(1) + 0;
t67 = t43 * t74;
t66 = t44 * pkin(1) + pkin(7) * t73 + 0;
t65 = t71 * pkin(7) + t68;
t64 = cos(t69);
t63 = sin(t70);
t62 = cos(t70) / 0.2e1;
t61 = sin(t69) / 0.2e1;
t60 = t42 * pkin(1) - pkin(7) * t72 + 0;
t47 = sin(qJ(2));
t56 = t62 + t64 / 0.2e1;
t23 = t42 * t56 + t44 * t47;
t31 = t61 - t63 / 0.2e1;
t49 = cos(qJ(2));
t24 = -t42 * t31 + t44 * t49;
t59 = t24 * pkin(2) + t23 * pkin(8) + t66;
t30 = t61 + t63 / 0.2e1;
t32 = t62 - t64 / 0.2e1;
t58 = t32 * pkin(2) - t30 * pkin(8) + t65;
t21 = t42 * t47 - t44 * t56;
t22 = t44 * t31 + t42 * t49;
t57 = t22 * pkin(2) + t21 * pkin(8) + t60;
t46 = sin(qJ(3));
t13 = t24 * t46 - t42 * t67;
t14 = t24 * t74 + t46 * t73;
t55 = t14 * pkin(3) + t13 * qJ(4) + t59;
t25 = t32 * t46 - t71 * t74;
t26 = t32 * t74 + t71 * t46;
t54 = t26 * pkin(3) + t25 * qJ(4) + t58;
t11 = t22 * t46 + t44 * t67;
t12 = t22 * t74 - t46 * t72;
t53 = t12 * pkin(3) + t11 * qJ(4) + t57;
t52 = t23 * pkin(4) + t14 * pkin(9) + t55;
t51 = -t30 * pkin(4) + t26 * pkin(9) + t54;
t50 = t21 * pkin(4) + t12 * pkin(9) + t53;
t48 = cos(qJ(5));
t45 = sin(qJ(5));
t6 = t25 * t45 - t30 * t48;
t5 = -t25 * t48 - t30 * t45;
t4 = t13 * t45 + t23 * t48;
t3 = -t13 * t48 + t23 * t45;
t2 = t11 * t45 + t21 * t48;
t1 = -t11 * t48 + t21 * t45;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t44, -t42, 0, 0; t42, t44, 0, 0; 0, 0, 1, t68; 0, 0, 0, 1; t24, -t23, t73, t66; t22, -t21, -t72, t60; t32, t30, t71, t65; 0, 0, 0, 1; t14, -t13, t23, t59; t12, -t11, t21, t57; t26, -t25, -t30, t58; 0, 0, 0, 1; t23, -t14, t13, t55; t21, -t12, t11, t53; -t30, -t26, t25, t54; 0, 0, 0, 1; t4, -t3, t14, t52; t2, -t1, t12, t50; t6, -t5, t26, t51; 0, 0, 0, 1; t4, t14, t3, t4 * pkin(5) + t3 * qJ(6) + t52; t2, t12, t1, t2 * pkin(5) + t1 * qJ(6) + t50; t6, t26, t5, t6 * pkin(5) + t5 * qJ(6) + t51; 0, 0, 0, 1;];
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
