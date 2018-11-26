% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:35
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRP11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:34:44
% EndTime: 2018-11-23 18:34:44
% DurationCPUTime: 0.29s
% Computational Cost: add. (1644->87), mult. (1692->105), div. (0->0), fcn. (1753->22), ass. (0->79)
t55 = pkin(6) - qJ(2);
t47 = cos(t55) / 0.2e1;
t54 = pkin(6) + qJ(2);
t50 = cos(t54);
t41 = t47 + t50 / 0.2e1;
t64 = sin(qJ(2));
t65 = sin(qJ(1));
t69 = cos(qJ(1));
t32 = -t65 * t41 - t69 * t64;
t56 = sin(pkin(7));
t58 = cos(pkin(7));
t57 = sin(pkin(6));
t98 = t65 * t57;
t83 = -t32 * t56 + t58 * t98;
t46 = sin(t54) / 0.2e1;
t49 = sin(t55);
t38 = t46 + t49 / 0.2e1;
t59 = cos(pkin(6));
t85 = -t38 * t56 + t59 * t58;
t102 = cos(qJ(4));
t97 = t69 * t57;
t96 = pkin(7) - qJ(3);
t95 = pkin(7) + qJ(3);
t94 = pkin(8) + 0;
t61 = sin(qJ(5));
t92 = pkin(5) * t61 + pkin(11);
t91 = t59 * pkin(9) + t94;
t90 = t69 * pkin(1) + pkin(9) * t98 + 0;
t89 = cos(t95);
t88 = sin(t96);
t87 = cos(t96) / 0.2e1;
t86 = sin(t95) / 0.2e1;
t30 = t69 * t41 - t65 * t64;
t84 = -t30 * t56 - t58 * t97;
t82 = t65 * pkin(1) - pkin(9) * t97 + 0;
t42 = t47 - t50 / 0.2e1;
t81 = t42 * pkin(2) + t85 * pkin(10) + t91;
t39 = t46 - t49 / 0.2e1;
t68 = cos(qJ(2));
t33 = -t65 * t39 + t69 * t68;
t80 = t33 * pkin(2) + t83 * pkin(10) + t90;
t37 = t86 - t88 / 0.2e1;
t40 = t87 - t89 / 0.2e1;
t67 = cos(qJ(3));
t21 = t38 * t37 + t59 * t40 + t42 * t67;
t79 = t21 * pkin(3) + t81;
t18 = t32 * t37 + t33 * t67 + t40 * t98;
t78 = t18 * pkin(3) + t80;
t77 = t87 + t89 / 0.2e1;
t76 = t86 + t88 / 0.2e1;
t75 = t57 * t76;
t63 = sin(qJ(3));
t20 = -t38 * t77 + t42 * t63 - t59 * t76;
t74 = t20 * pkin(11) + t79;
t17 = -t32 * t77 + t33 * t63 - t65 * t75;
t73 = t17 * pkin(11) + t78;
t31 = t69 * t39 + t65 * t68;
t72 = t31 * pkin(2) + t84 * pkin(10) + t82;
t16 = t30 * t37 + t31 * t67 - t40 * t97;
t71 = t16 * pkin(3) + t72;
t15 = -t30 * t77 + t31 * t63 + t69 * t75;
t70 = t15 * pkin(11) + t71;
t66 = cos(qJ(5));
t62 = sin(qJ(4));
t60 = -qJ(6) - pkin(12);
t48 = t66 * pkin(5) + pkin(4);
t12 = t21 * t102 + t85 * t62;
t11 = -t85 * t102 + t21 * t62;
t10 = t18 * t102 + t83 * t62;
t9 = -t83 * t102 + t18 * t62;
t8 = t16 * t102 + t84 * t62;
t7 = -t84 * t102 + t16 * t62;
t6 = t12 * t66 + t20 * t61;
t5 = -t12 * t61 + t20 * t66;
t4 = t10 * t66 + t17 * t61;
t3 = -t10 * t61 + t17 * t66;
t2 = t15 * t61 + t8 * t66;
t1 = t15 * t66 - t8 * t61;
t13 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t69, -t65, 0, 0; t65, t69, 0, 0; 0, 0, 1, t94; 0, 0, 0, 1; t33, t32, t98, t90; t31, t30, -t97, t82; t42, t38, t59, t91; 0, 0, 0, 1; t18, -t17, t83, t80; t16, -t15, t84, t72; t21, -t20, t85, t81; 0, 0, 0, 1; t10, -t9, t17, t73; t8, -t7, t15, t70; t12, -t11, t20, t74; 0, 0, 0, 1; t4, t3, t9, t10 * pkin(4) + t9 * pkin(12) + t73; t2, t1, t7, t8 * pkin(4) + t7 * pkin(12) + t70; t6, t5, t11, t12 * pkin(4) + t11 * pkin(12) + t74; 0, 0, 0, 1; t4, t3, t9, t10 * t48 + t92 * t17 - t9 * t60 + t78; t2, t1, t7, t92 * t15 + t8 * t48 - t7 * t60 + t71; t6, t5, t11, -t11 * t60 + t12 * t48 + t92 * t20 + t79; 0, 0, 0, 1;];
T_ges = t13;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
