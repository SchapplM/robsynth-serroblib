% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRRP11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:30:27
% EndTime: 2018-11-23 16:30:27
% DurationCPUTime: 0.30s
% Computational Cost: add. (1644->87), mult. (1692->106), div. (0->0), fcn. (1753->22), ass. (0->80)
t55 = pkin(6) - pkin(12);
t47 = cos(t55) / 0.2e1;
t54 = pkin(6) + pkin(12);
t50 = cos(t54);
t39 = t47 + t50 / 0.2e1;
t56 = sin(pkin(12));
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t32 = -t66 * t39 - t69 * t56;
t57 = sin(pkin(7));
t60 = cos(pkin(7));
t58 = sin(pkin(6));
t99 = t66 * t58;
t83 = -t32 * t57 + t60 * t99;
t46 = sin(t54) / 0.2e1;
t49 = sin(t55);
t37 = t46 + t49 / 0.2e1;
t61 = cos(pkin(6));
t85 = -t37 * t57 + t61 * t60;
t103 = cos(qJ(4));
t98 = t69 * t58;
t97 = qJ(2) * t58;
t96 = pkin(7) - qJ(3);
t95 = pkin(7) + qJ(3);
t94 = pkin(8) + 0;
t63 = sin(qJ(5));
t92 = pkin(5) * t63 + pkin(10);
t91 = t61 * qJ(2) + t94;
t90 = t69 * pkin(1) + t66 * t97 + 0;
t89 = cos(t95);
t88 = sin(t96);
t87 = cos(t96) / 0.2e1;
t86 = sin(t95) / 0.2e1;
t30 = t69 * t39 - t66 * t56;
t84 = -t30 * t57 - t60 * t98;
t82 = t66 * pkin(1) - t69 * t97 + 0;
t40 = t47 - t50 / 0.2e1;
t81 = t40 * pkin(2) + t85 * pkin(9) + t91;
t38 = t46 - t49 / 0.2e1;
t59 = cos(pkin(12));
t33 = -t66 * t38 + t69 * t59;
t80 = t33 * pkin(2) + t83 * pkin(9) + t90;
t41 = t86 - t88 / 0.2e1;
t42 = t87 - t89 / 0.2e1;
t68 = cos(qJ(3));
t21 = t37 * t41 + t40 * t68 + t61 * t42;
t79 = t21 * pkin(3) + t81;
t18 = t32 * t41 + t33 * t68 + t42 * t99;
t78 = t18 * pkin(3) + t80;
t77 = t87 + t89 / 0.2e1;
t76 = t86 + t88 / 0.2e1;
t75 = t58 * t76;
t65 = sin(qJ(3));
t20 = -t37 * t77 + t40 * t65 - t61 * t76;
t74 = t20 * pkin(10) + t79;
t17 = -t32 * t77 + t33 * t65 - t66 * t75;
t73 = t17 * pkin(10) + t78;
t31 = t69 * t38 + t66 * t59;
t72 = t31 * pkin(2) + t84 * pkin(9) + t82;
t16 = t30 * t41 + t31 * t68 - t42 * t98;
t71 = t16 * pkin(3) + t72;
t15 = -t30 * t77 + t31 * t65 + t69 * t75;
t70 = t15 * pkin(10) + t71;
t67 = cos(qJ(5));
t64 = sin(qJ(4));
t62 = -qJ(6) - pkin(11);
t48 = t67 * pkin(5) + pkin(4);
t12 = t21 * t103 + t85 * t64;
t11 = -t85 * t103 + t21 * t64;
t10 = t18 * t103 + t83 * t64;
t9 = -t83 * t103 + t18 * t64;
t8 = t16 * t103 + t84 * t64;
t7 = -t84 * t103 + t16 * t64;
t6 = t12 * t67 + t20 * t63;
t5 = -t12 * t63 + t20 * t67;
t4 = t10 * t67 + t17 * t63;
t3 = -t10 * t63 + t17 * t67;
t2 = t15 * t63 + t8 * t67;
t1 = t15 * t67 - t8 * t63;
t13 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t69, -t66, 0, 0; t66, t69, 0, 0; 0, 0, 1, t94; 0, 0, 0, 1; t33, t32, t99, t90; t31, t30, -t98, t82; t40, t37, t61, t91; 0, 0, 0, 1; t18, -t17, t83, t80; t16, -t15, t84, t72; t21, -t20, t85, t81; 0, 0, 0, 1; t10, -t9, t17, t73; t8, -t7, t15, t70; t12, -t11, t20, t74; 0, 0, 0, 1; t4, t3, t9, t10 * pkin(4) + t9 * pkin(11) + t73; t2, t1, t7, t8 * pkin(4) + t7 * pkin(11) + t70; t6, t5, t11, t12 * pkin(4) + t11 * pkin(11) + t74; 0, 0, 0, 1; t4, t3, t9, t10 * t48 + t92 * t17 - t9 * t62 + t78; t2, t1, t7, t92 * t15 + t8 * t48 - t7 * t62 + t71; t6, t5, t11, -t11 * t62 + t12 * t48 + t92 * t20 + t79; 0, 0, 0, 1;];
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
