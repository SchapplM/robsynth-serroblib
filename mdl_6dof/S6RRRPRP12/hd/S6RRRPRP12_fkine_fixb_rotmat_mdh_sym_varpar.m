% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2018-11-23 17:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPRP12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:49:43
% EndTime: 2018-11-23 17:49:43
% DurationCPUTime: 0.15s
% Computational Cost: add. (611->67), mult. (693->70), div. (0->0), fcn. (780->14), ass. (0->53)
t74 = cos(qJ(3));
t42 = sin(pkin(6));
t46 = sin(qJ(1));
t73 = t46 * t42;
t49 = cos(qJ(1));
t72 = t49 * t42;
t71 = cos(pkin(6));
t70 = pkin(6) - qJ(2);
t69 = pkin(6) + qJ(2);
t68 = pkin(7) + 0;
t67 = t42 * t74;
t66 = t71 * pkin(8) + t68;
t65 = t49 * pkin(1) + pkin(8) * t73 + 0;
t64 = cos(t69);
t63 = sin(t70);
t62 = cos(t70) / 0.2e1;
t61 = sin(t69) / 0.2e1;
t60 = t46 * pkin(1) - pkin(8) * t72 + 0;
t30 = t61 + t63 / 0.2e1;
t32 = t62 - t64 / 0.2e1;
t59 = t32 * pkin(2) - t30 * pkin(9) + t66;
t45 = sin(qJ(2));
t56 = t62 + t64 / 0.2e1;
t25 = t49 * t45 + t46 * t56;
t31 = t61 - t63 / 0.2e1;
t48 = cos(qJ(2));
t26 = -t46 * t31 + t49 * t48;
t58 = t26 * pkin(2) + t25 * pkin(9) + t65;
t23 = t46 * t45 - t49 * t56;
t24 = t49 * t31 + t46 * t48;
t57 = t24 * pkin(2) + t23 * pkin(9) + t60;
t44 = sin(qJ(3));
t21 = t32 * t44 - t71 * t74;
t22 = t32 * t74 + t71 * t44;
t55 = t22 * pkin(3) + t21 * qJ(4) + t59;
t13 = t26 * t44 - t46 * t67;
t14 = t26 * t74 + t44 * t73;
t54 = t14 * pkin(3) + t13 * qJ(4) + t58;
t11 = t24 * t44 + t49 * t67;
t12 = t24 * t74 - t44 * t72;
t53 = t12 * pkin(3) + t11 * qJ(4) + t57;
t52 = -t30 * pkin(4) + t22 * pkin(10) + t55;
t51 = t25 * pkin(4) + t14 * pkin(10) + t54;
t50 = t23 * pkin(4) + t12 * pkin(10) + t53;
t47 = cos(qJ(5));
t43 = sin(qJ(5));
t6 = t21 * t43 - t30 * t47;
t5 = -t21 * t47 - t30 * t43;
t4 = t13 * t43 + t25 * t47;
t3 = -t13 * t47 + t25 * t43;
t2 = t11 * t43 + t23 * t47;
t1 = -t11 * t47 + t23 * t43;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t49, -t46, 0, 0; t46, t49, 0, 0; 0, 0, 1, t68; 0, 0, 0, 1; t26, -t25, t73, t65; t24, -t23, -t72, t60; t32, t30, t71, t66; 0, 0, 0, 1; t14, -t13, t25, t58; t12, -t11, t23, t57; t22, -t21, -t30, t59; 0, 0, 0, 1; t25, -t14, t13, t54; t23, -t12, t11, t53; -t30, -t22, t21, t55; 0, 0, 0, 1; t4, -t3, t14, t51; t2, -t1, t12, t50; t6, -t5, t22, t52; 0, 0, 0, 1; t4, t14, t3, t4 * pkin(5) + t3 * qJ(6) + t51; t2, t12, t1, t2 * pkin(5) + t1 * qJ(6) + t50; t6, t22, t5, t6 * pkin(5) + t5 * qJ(6) + t52; 0, 0, 0, 1;];
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
