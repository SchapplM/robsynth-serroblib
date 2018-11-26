% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRRPR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:25:49
% EndTime: 2018-11-23 15:25:50
% DurationCPUTime: 0.32s
% Computational Cost: add. (1656->94), mult. (1692->117), div. (0->0), fcn. (1753->24), ass. (0->76)
t52 = pkin(6) - qJ(2);
t41 = cos(t52) / 0.2e1;
t51 = pkin(6) + qJ(2);
t46 = cos(t51);
t35 = t41 + t46 / 0.2e1;
t54 = sin(pkin(12));
t58 = cos(pkin(12));
t64 = sin(qJ(2));
t26 = -t54 * t35 - t58 * t64;
t55 = sin(pkin(7));
t59 = cos(pkin(7));
t56 = sin(pkin(6));
t96 = t54 * t56;
t80 = -t26 * t55 + t59 * t96;
t40 = sin(t51) / 0.2e1;
t44 = sin(t52);
t32 = t40 + t44 / 0.2e1;
t60 = cos(pkin(6));
t82 = -t32 * t55 + t60 * t59;
t99 = cos(qJ(4));
t95 = t58 * t56;
t93 = pkin(7) - qJ(3);
t92 = pkin(7) + qJ(3);
t90 = qJ(1) + 0;
t53 = sin(pkin(13));
t89 = pkin(5) * t53 + pkin(10);
t88 = t58 * pkin(1) + pkin(8) * t96 + 0;
t87 = t60 * pkin(8) + t90;
t86 = cos(t92);
t85 = sin(t93);
t84 = cos(t93) / 0.2e1;
t83 = sin(t92) / 0.2e1;
t24 = t58 * t35 - t54 * t64;
t81 = -t24 * t55 - t59 * t95;
t79 = t54 * pkin(1) - pkin(8) * t95 + 0;
t33 = t40 - t44 / 0.2e1;
t66 = cos(qJ(2));
t27 = -t54 * t33 + t58 * t66;
t78 = t27 * pkin(2) + pkin(9) * t80 + t88;
t36 = t41 - t46 / 0.2e1;
t77 = t36 * pkin(2) + pkin(9) * t82 + t87;
t31 = t83 - t85 / 0.2e1;
t34 = t84 - t86 / 0.2e1;
t65 = cos(qJ(3));
t12 = t26 * t31 + t27 * t65 + t34 * t96;
t76 = t12 * pkin(3) + t78;
t15 = t32 * t31 + t60 * t34 + t36 * t65;
t75 = t15 * pkin(3) + t77;
t74 = t84 + t86 / 0.2e1;
t73 = t83 + t85 / 0.2e1;
t72 = t56 * t73;
t63 = sin(qJ(3));
t11 = -t26 * t74 + t27 * t63 - t54 * t72;
t71 = t11 * pkin(10) + t76;
t14 = -t32 * t74 + t36 * t63 - t60 * t73;
t70 = t14 * pkin(10) + t75;
t25 = t58 * t33 + t54 * t66;
t69 = t25 * pkin(2) + t81 * pkin(9) + t79;
t10 = t24 * t31 + t25 * t65 - t34 * t95;
t68 = t10 * pkin(3) + t69;
t9 = -t24 * t74 + t25 * t63 + t58 * t72;
t67 = t9 * pkin(10) + t68;
t62 = sin(qJ(4));
t61 = -pkin(11) - qJ(5);
t57 = cos(pkin(13));
t50 = pkin(13) + qJ(6);
t45 = cos(t50);
t43 = sin(t50);
t42 = t57 * pkin(5) + pkin(4);
t6 = t15 * t99 + t82 * t62;
t5 = t15 * t62 - t82 * t99;
t4 = t12 * t99 + t80 * t62;
t3 = t12 * t62 - t80 * t99;
t2 = t10 * t99 + t81 * t62;
t1 = t10 * t62 - t81 * t99;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t58, -t54, 0, 0; t54, t58, 0, 0; 0, 0, 1, t90; 0, 0, 0, 1; t27, t26, t96, t88; t25, t24, -t95, t79; t36, t32, t60, t87; 0, 0, 0, 1; t12, -t11, t80, t78; t10, -t9, t81, t69; t15, -t14, t82, t77; 0, 0, 0, 1; t4, -t3, t11, t71; t2, -t1, t9, t67; t6, -t5, t14, t70; 0, 0, 0, 1; t11 * t53 + t4 * t57, t11 * t57 - t4 * t53, t3, t4 * pkin(4) + t3 * qJ(5) + t71; t2 * t57 + t9 * t53, -t2 * t53 + t9 * t57, t1, t2 * pkin(4) + t1 * qJ(5) + t67; t14 * t53 + t6 * t57, t14 * t57 - t6 * t53, t5, t6 * pkin(4) + t5 * qJ(5) + t70; 0, 0, 0, 1; t11 * t43 + t4 * t45, t11 * t45 - t4 * t43, t3, t89 * t11 - t3 * t61 + t4 * t42 + t76; t2 * t45 + t9 * t43, -t2 * t43 + t9 * t45, t1, -t1 * t61 + t2 * t42 + t89 * t9 + t68; t14 * t43 + t6 * t45, t14 * t45 - t6 * t43, t5, t89 * t14 + t6 * t42 - t5 * t61 + t75; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
