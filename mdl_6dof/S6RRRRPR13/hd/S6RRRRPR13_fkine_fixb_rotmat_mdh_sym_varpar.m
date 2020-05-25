% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRRPR13_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:22:37
% EndTime: 2018-11-23 18:22:37
% DurationCPUTime: 0.19s
% Computational Cost: add. (752->76), mult. (877->82), div. (0->0), fcn. (997->16), ass. (0->59)
t68 = pkin(6) + qJ(2);
t79 = sin(t68) / 0.2e1;
t78 = pkin(10) - pkin(11);
t69 = pkin(6) - qJ(2);
t62 = sin(t69);
t31 = t79 - t62 / 0.2e1;
t46 = sin(qJ(1));
t48 = cos(qJ(2));
t49 = cos(qJ(1));
t24 = t49 * t31 + t46 * t48;
t44 = sin(qJ(3));
t41 = sin(pkin(6));
t74 = cos(qJ(3));
t66 = t41 * t74;
t12 = t24 * t44 + t49 * t66;
t77 = t12 * pkin(10);
t26 = -t46 * t31 + t49 * t48;
t14 = t26 * t44 - t46 * t66;
t76 = t14 * pkin(10);
t61 = cos(t69) / 0.2e1;
t63 = cos(t68);
t32 = t61 - t63 / 0.2e1;
t70 = cos(pkin(6));
t21 = t32 * t44 - t70 * t74;
t75 = t21 * pkin(10);
t73 = cos(qJ(4));
t72 = t46 * t41;
t71 = t49 * t41;
t67 = pkin(7) + 0;
t65 = t70 * pkin(8) + t67;
t64 = t49 * pkin(1) + pkin(8) * t72 + 0;
t60 = t46 * pkin(1) - pkin(8) * t71 + 0;
t30 = t79 + t62 / 0.2e1;
t59 = t32 * pkin(2) - t30 * pkin(9) + t65;
t45 = sin(qJ(2));
t54 = t61 + t63 / 0.2e1;
t25 = t49 * t45 + t46 * t54;
t58 = t26 * pkin(2) + t25 * pkin(9) + t64;
t22 = t32 * t74 + t70 * t44;
t57 = t22 * pkin(3) + t59;
t15 = t26 * t74 + t44 * t72;
t56 = t15 * pkin(3) + t58;
t23 = t46 * t45 - t49 * t54;
t55 = t24 * pkin(2) + t23 * pkin(9) + t60;
t13 = t24 * t74 - t44 * t71;
t53 = t13 * pkin(3) + t55;
t43 = sin(qJ(4));
t8 = t22 * t43 + t30 * t73;
t9 = t22 * t73 - t30 * t43;
t52 = t9 * pkin(4) + t8 * qJ(5) + t57;
t5 = t15 * t43 - t25 * t73;
t6 = t15 * t73 + t25 * t43;
t51 = t6 * pkin(4) + t5 * qJ(5) + t56;
t3 = t13 * t43 - t23 * t73;
t4 = t13 * t73 + t23 * t43;
t50 = t4 * pkin(4) + t3 * qJ(5) + t53;
t47 = cos(qJ(6));
t42 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t49, -t46, 0, 0; t46, t49, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t26, -t25, t72, t64; t24, -t23, -t71, t60; t32, t30, t70, t65; 0, 0, 0, 1; t15, -t14, t25, t58; t13, -t12, t23, t55; t22, -t21, -t30, t59; 0, 0, 0, 1; t6, -t5, t14, t56 + t76; t4, -t3, t12, t53 + t77; t9, -t8, t21, t57 + t75; 0, 0, 0, 1; t6, t14, t5, t51 + t76; t4, t12, t3, t50 + t77; t9, t21, t8, t52 + t75; 0, 0, 0, 1; t5 * t42 + t6 * t47, -t6 * t42 + t5 * t47, -t14, t6 * pkin(5) + t78 * t14 + t51; t3 * t42 + t4 * t47, t3 * t47 - t4 * t42, -t12, t4 * pkin(5) + t78 * t12 + t50; t8 * t42 + t9 * t47, -t9 * t42 + t8 * t47, -t21, t9 * pkin(5) + t78 * t21 + t52; 0, 0, 0, 1;];
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
