% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRRPRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:19:13
% EndTime: 2018-11-23 15:19:13
% DurationCPUTime: 0.29s
% Computational Cost: add. (1408->84), mult. (1439->101), div. (0->0), fcn. (1483->22), ass. (0->71)
t54 = pkin(6) - qJ(2);
t47 = cos(t54) / 0.2e1;
t53 = pkin(6) + qJ(2);
t49 = cos(t53);
t41 = t47 + t49 / 0.2e1;
t55 = sin(pkin(12));
t58 = cos(pkin(12));
t64 = sin(qJ(2));
t32 = -t41 * t55 - t58 * t64;
t56 = sin(pkin(7));
t59 = cos(pkin(7));
t57 = sin(pkin(6));
t96 = t55 * t57;
t22 = -t32 * t56 + t59 * t96;
t46 = sin(t53) / 0.2e1;
t48 = sin(t54);
t39 = t46 + t48 / 0.2e1;
t60 = cos(pkin(6));
t29 = -t39 * t56 + t59 * t60;
t95 = t58 * t57;
t93 = pkin(7) - qJ(3);
t92 = pkin(7) + qJ(3);
t90 = qJ(1) + 0;
t89 = pkin(1) * t58 + pkin(8) * t96 + 0;
t88 = pkin(8) * t60 + t90;
t87 = cos(t92);
t86 = sin(t93);
t85 = cos(t93) / 0.2e1;
t84 = sin(t92) / 0.2e1;
t30 = t41 * t58 - t55 * t64;
t21 = -t30 * t56 - t59 * t95;
t83 = pkin(1) * t55 - pkin(8) * t95 + 0;
t40 = t46 - t48 / 0.2e1;
t68 = cos(qJ(2));
t33 = -t40 * t55 + t58 * t68;
t82 = t33 * pkin(2) + pkin(9) * t22 + t89;
t42 = t47 - t49 / 0.2e1;
t81 = t42 * pkin(2) + pkin(9) * t29 + t88;
t80 = t85 - t87 / 0.2e1;
t79 = t85 + t87 / 0.2e1;
t78 = t84 + t86 / 0.2e1;
t77 = t57 * t80;
t76 = t57 * t78;
t63 = sin(qJ(3));
t13 = -t32 * t79 + t33 * t63 - t55 * t76;
t38 = t84 - t86 / 0.2e1;
t67 = cos(qJ(3));
t14 = t32 * t38 + t33 * t67 + t55 * t77;
t75 = pkin(3) * t14 + qJ(4) * t13 + t82;
t17 = -t39 * t79 + t42 * t63 - t60 * t78;
t18 = t38 * t39 + t42 * t67 + t60 * t80;
t74 = pkin(3) * t18 + qJ(4) * t17 + t81;
t31 = t40 * t58 + t55 * t68;
t73 = pkin(2) * t31 + pkin(9) * t21 + t83;
t72 = pkin(4) * t22 + pkin(10) * t14 + t75;
t71 = pkin(4) * t29 + pkin(10) * t18 + t74;
t11 = -t30 * t79 + t31 * t63 + t58 * t76;
t12 = t30 * t38 + t31 * t67 - t58 * t77;
t70 = pkin(3) * t12 + t11 * qJ(4) + t73;
t69 = pkin(4) * t21 + t12 * pkin(10) + t70;
t66 = cos(qJ(5));
t65 = cos(qJ(6));
t62 = sin(qJ(5));
t61 = sin(qJ(6));
t6 = t17 * t62 + t29 * t66;
t5 = -t17 * t66 + t29 * t62;
t4 = t13 * t62 + t22 * t66;
t3 = -t13 * t66 + t22 * t62;
t2 = t11 * t62 + t21 * t66;
t1 = -t11 * t66 + t21 * t62;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t58, -t55, 0, 0; t55, t58, 0, 0; 0, 0, 1, t90; 0, 0, 0, 1; t33, t32, t96, t89; t31, t30, -t95, t83; t42, t39, t60, t88; 0, 0, 0, 1; t14, -t13, t22, t82; t12, -t11, t21, t73; t18, -t17, t29, t81; 0, 0, 0, 1; t22, -t14, t13, t75; t21, -t12, t11, t70; t29, -t18, t17, t74; 0, 0, 0, 1; t4, -t3, t14, t72; t2, -t1, t12, t69; t6, -t5, t18, t71; 0, 0, 0, 1; t14 * t61 + t4 * t65, t14 * t65 - t4 * t61, t3, pkin(5) * t4 + pkin(11) * t3 + t72; t12 * t61 + t2 * t65, t12 * t65 - t2 * t61, t1, t2 * pkin(5) + t1 * pkin(11) + t69; t18 * t61 + t6 * t65, t18 * t65 - t6 * t61, t5, pkin(5) * t6 + pkin(11) * t5 + t71; 0, 0, 0, 1;];
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
