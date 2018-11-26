% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2018-11-23 14:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PPRRRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:50:43
% EndTime: 2018-11-23 14:50:43
% DurationCPUTime: 0.31s
% Computational Cost: add. (1644->87), mult. (1692->106), div. (0->0), fcn. (1753->22), ass. (0->80)
t57 = sin(pkin(11));
t59 = sin(pkin(6));
t100 = t57 * t59;
t55 = pkin(6) - pkin(12);
t47 = cos(t55) / 0.2e1;
t54 = pkin(6) + pkin(12);
t50 = cos(t54);
t39 = t47 + t50 / 0.2e1;
t56 = sin(pkin(12));
t61 = cos(pkin(11));
t32 = -t57 * t39 - t61 * t56;
t58 = sin(pkin(7));
t62 = cos(pkin(7));
t83 = t62 * t100 - t32 * t58;
t46 = sin(t54) / 0.2e1;
t49 = sin(t55);
t37 = t46 + t49 / 0.2e1;
t63 = cos(pkin(6));
t85 = -t37 * t58 + t63 * t62;
t103 = cos(qJ(4));
t99 = t61 * t59;
t97 = qJ(2) * t59;
t96 = pkin(7) - qJ(3);
t95 = pkin(7) + qJ(3);
t93 = qJ(1) + 0;
t65 = sin(qJ(5));
t92 = pkin(5) * t65 + pkin(9);
t91 = t61 * pkin(1) + t57 * t97 + 0;
t90 = t63 * qJ(2) + t93;
t89 = cos(t95);
t88 = sin(t96);
t87 = cos(t96) / 0.2e1;
t86 = sin(t95) / 0.2e1;
t30 = t61 * t39 - t57 * t56;
t84 = -t30 * t58 - t62 * t99;
t82 = t57 * pkin(1) - t61 * t97 + 0;
t38 = t46 - t49 / 0.2e1;
t60 = cos(pkin(12));
t33 = -t57 * t38 + t61 * t60;
t81 = t33 * pkin(2) + t83 * pkin(8) + t91;
t40 = t47 - t50 / 0.2e1;
t80 = t40 * pkin(2) + t85 * pkin(8) + t90;
t41 = t86 - t88 / 0.2e1;
t42 = t87 - t89 / 0.2e1;
t69 = cos(qJ(3));
t18 = t42 * t100 + t32 * t41 + t33 * t69;
t79 = t18 * pkin(3) + t81;
t21 = t37 * t41 + t40 * t69 + t63 * t42;
t78 = t21 * pkin(3) + t80;
t77 = t87 + t89 / 0.2e1;
t76 = t86 + t88 / 0.2e1;
t75 = t59 * t76;
t67 = sin(qJ(3));
t17 = -t32 * t77 + t33 * t67 - t57 * t75;
t74 = t17 * pkin(9) + t79;
t20 = -t37 * t77 + t40 * t67 - t63 * t76;
t73 = t20 * pkin(9) + t78;
t31 = t61 * t38 + t57 * t60;
t72 = t31 * pkin(2) + t84 * pkin(8) + t82;
t16 = t30 * t41 + t31 * t69 - t42 * t99;
t71 = t16 * pkin(3) + t72;
t15 = -t30 * t77 + t31 * t67 + t61 * t75;
t70 = t15 * pkin(9) + t71;
t68 = cos(qJ(5));
t66 = sin(qJ(4));
t64 = -qJ(6) - pkin(10);
t48 = t68 * pkin(5) + pkin(4);
t12 = t21 * t103 + t85 * t66;
t11 = -t85 * t103 + t21 * t66;
t10 = t18 * t103 + t83 * t66;
t9 = -t83 * t103 + t18 * t66;
t8 = t16 * t103 + t84 * t66;
t7 = -t84 * t103 + t16 * t66;
t6 = t12 * t68 + t20 * t65;
t5 = -t12 * t65 + t20 * t68;
t4 = t10 * t68 + t17 * t65;
t3 = -t10 * t65 + t17 * t68;
t2 = t15 * t65 + t8 * t68;
t1 = t15 * t68 - t8 * t65;
t13 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t61, -t57, 0, 0; t57, t61, 0, 0; 0, 0, 1, t93; 0, 0, 0, 1; t33, t32, t100, t91; t31, t30, -t99, t82; t40, t37, t63, t90; 0, 0, 0, 1; t18, -t17, t83, t81; t16, -t15, t84, t72; t21, -t20, t85, t80; 0, 0, 0, 1; t10, -t9, t17, t74; t8, -t7, t15, t70; t12, -t11, t20, t73; 0, 0, 0, 1; t4, t3, t9, t10 * pkin(4) + t9 * pkin(10) + t74; t2, t1, t7, t8 * pkin(4) + t7 * pkin(10) + t70; t6, t5, t11, t12 * pkin(4) + t11 * pkin(10) + t73; 0, 0, 0, 1; t4, t3, t9, t10 * t48 + t92 * t17 - t9 * t64 + t79; t2, t1, t7, t92 * t15 + t8 * t48 - t7 * t64 + t71; t6, t5, t11, -t11 * t64 + t12 * t48 + t92 * t20 + t78; 0, 0, 0, 1;];
T_ges = t13;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
