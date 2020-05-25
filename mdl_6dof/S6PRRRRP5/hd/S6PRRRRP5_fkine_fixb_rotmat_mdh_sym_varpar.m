% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRRRRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:30:19
% EndTime: 2018-11-23 15:30:19
% DurationCPUTime: 0.31s
% Computational Cost: add. (1644->87), mult. (1692->105), div. (0->0), fcn. (1753->22), ass. (0->79)
t55 = pkin(6) - qJ(2);
t47 = cos(t55) / 0.2e1;
t54 = pkin(6) + qJ(2);
t50 = cos(t54);
t41 = t47 + t50 / 0.2e1;
t56 = sin(pkin(12));
t59 = cos(pkin(12));
t66 = sin(qJ(2));
t32 = -t56 * t41 - t59 * t66;
t57 = sin(pkin(7));
t60 = cos(pkin(7));
t58 = sin(pkin(6));
t99 = t56 * t58;
t83 = -t32 * t57 + t60 * t99;
t46 = sin(t54) / 0.2e1;
t49 = sin(t55);
t38 = t46 + t49 / 0.2e1;
t61 = cos(pkin(6));
t85 = -t38 * t57 + t61 * t60;
t102 = cos(qJ(4));
t98 = t59 * t58;
t96 = pkin(7) - qJ(3);
t95 = pkin(7) + qJ(3);
t93 = qJ(1) + 0;
t63 = sin(qJ(5));
t92 = pkin(5) * t63 + pkin(10);
t91 = t59 * pkin(1) + pkin(8) * t99 + 0;
t90 = t61 * pkin(8) + t93;
t89 = cos(t95);
t88 = sin(t96);
t87 = cos(t96) / 0.2e1;
t86 = sin(t95) / 0.2e1;
t30 = t59 * t41 - t56 * t66;
t84 = -t30 * t57 - t60 * t98;
t82 = t56 * pkin(1) - pkin(8) * t98 + 0;
t39 = t46 - t49 / 0.2e1;
t69 = cos(qJ(2));
t33 = -t56 * t39 + t59 * t69;
t81 = t33 * pkin(2) + t83 * pkin(9) + t91;
t42 = t47 - t50 / 0.2e1;
t80 = t42 * pkin(2) + t85 * pkin(9) + t90;
t37 = t86 - t88 / 0.2e1;
t40 = t87 - t89 / 0.2e1;
t68 = cos(qJ(3));
t18 = t32 * t37 + t33 * t68 + t40 * t99;
t79 = t18 * pkin(3) + t81;
t21 = t38 * t37 + t61 * t40 + t42 * t68;
t78 = t21 * pkin(3) + t80;
t77 = t87 + t89 / 0.2e1;
t76 = t86 + t88 / 0.2e1;
t75 = t58 * t76;
t65 = sin(qJ(3));
t17 = -t32 * t77 + t33 * t65 - t56 * t75;
t74 = t17 * pkin(10) + t79;
t20 = -t38 * t77 + t42 * t65 - t61 * t76;
t73 = t20 * pkin(10) + t78;
t31 = t59 * t39 + t56 * t69;
t72 = t31 * pkin(2) + t84 * pkin(9) + t82;
t16 = t30 * t37 + t31 * t68 - t40 * t98;
t71 = t16 * pkin(3) + t72;
t15 = -t30 * t77 + t31 * t65 + t59 * t75;
t70 = t15 * pkin(10) + t71;
t67 = cos(qJ(5));
t64 = sin(qJ(4));
t62 = -qJ(6) - pkin(11);
t48 = t67 * pkin(5) + pkin(4);
t12 = t21 * t102 + t85 * t64;
t11 = -t85 * t102 + t21 * t64;
t10 = t18 * t102 + t83 * t64;
t9 = -t83 * t102 + t18 * t64;
t8 = t16 * t102 + t84 * t64;
t7 = -t84 * t102 + t16 * t64;
t6 = t12 * t67 + t20 * t63;
t5 = -t12 * t63 + t20 * t67;
t4 = t10 * t67 + t17 * t63;
t3 = -t10 * t63 + t17 * t67;
t2 = t15 * t63 + t8 * t67;
t1 = t15 * t67 - t8 * t63;
t13 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t59, -t56, 0, 0; t56, t59, 0, 0; 0, 0, 1, t93; 0, 0, 0, 1; t33, t32, t99, t91; t31, t30, -t98, t82; t42, t38, t61, t90; 0, 0, 0, 1; t18, -t17, t83, t81; t16, -t15, t84, t72; t21, -t20, t85, t80; 0, 0, 0, 1; t10, -t9, t17, t74; t8, -t7, t15, t70; t12, -t11, t20, t73; 0, 0, 0, 1; t4, t3, t9, t10 * pkin(4) + t9 * pkin(11) + t74; t2, t1, t7, t8 * pkin(4) + t7 * pkin(11) + t70; t6, t5, t11, t12 * pkin(4) + t11 * pkin(11) + t73; 0, 0, 0, 1; t4, t3, t9, t10 * t48 + t92 * t17 - t9 * t62 + t79; t2, t1, t7, t92 * t15 + t8 * t48 - t7 * t62 + t71; t6, t5, t11, -t11 * t62 + t12 * t48 + t92 * t20 + t78; 0, 0, 0, 1;];
T_ges = t13;
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
