% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRPR14_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:23:33
% EndTime: 2018-11-23 18:23:33
% DurationCPUTime: 0.30s
% Computational Cost: add. (1656->94), mult. (1692->117), div. (0->0), fcn. (1753->24), ass. (0->76)
t52 = pkin(6) - qJ(2);
t41 = cos(t52) / 0.2e1;
t51 = pkin(6) + qJ(2);
t46 = cos(t51);
t35 = t41 + t46 / 0.2e1;
t62 = sin(qJ(2));
t63 = sin(qJ(1));
t66 = cos(qJ(1));
t26 = -t63 * t35 - t66 * t62;
t54 = sin(pkin(7));
t57 = cos(pkin(7));
t55 = sin(pkin(6));
t95 = t63 * t55;
t80 = -t26 * t54 + t57 * t95;
t40 = sin(t51) / 0.2e1;
t44 = sin(t52);
t32 = t40 + t44 / 0.2e1;
t58 = cos(pkin(6));
t82 = -t32 * t54 + t58 * t57;
t99 = cos(qJ(4));
t94 = t66 * t55;
t93 = pkin(7) - qJ(3);
t92 = pkin(7) + qJ(3);
t91 = pkin(8) + 0;
t53 = sin(pkin(13));
t89 = pkin(5) * t53 + pkin(11);
t88 = t58 * pkin(9) + t91;
t87 = t66 * pkin(1) + pkin(9) * t95 + 0;
t86 = cos(t92);
t85 = sin(t93);
t84 = cos(t93) / 0.2e1;
t83 = sin(t92) / 0.2e1;
t24 = t66 * t35 - t63 * t62;
t81 = -t24 * t54 - t57 * t94;
t79 = t63 * pkin(1) - pkin(9) * t94 + 0;
t36 = t41 - t46 / 0.2e1;
t78 = t36 * pkin(2) + pkin(10) * t82 + t88;
t33 = t40 - t44 / 0.2e1;
t65 = cos(qJ(2));
t27 = -t63 * t33 + t66 * t65;
t77 = t27 * pkin(2) + pkin(10) * t80 + t87;
t31 = t83 - t85 / 0.2e1;
t34 = t84 - t86 / 0.2e1;
t64 = cos(qJ(3));
t15 = t32 * t31 + t58 * t34 + t36 * t64;
t76 = t15 * pkin(3) + t78;
t12 = t26 * t31 + t27 * t64 + t34 * t95;
t75 = t12 * pkin(3) + t77;
t74 = t84 + t86 / 0.2e1;
t73 = t83 + t85 / 0.2e1;
t72 = t55 * t73;
t61 = sin(qJ(3));
t14 = -t32 * t74 + t36 * t61 - t58 * t73;
t71 = t14 * pkin(11) + t76;
t11 = -t26 * t74 + t27 * t61 - t63 * t72;
t70 = t11 * pkin(11) + t75;
t25 = t66 * t33 + t63 * t65;
t69 = t25 * pkin(2) + t81 * pkin(10) + t79;
t10 = t24 * t31 + t25 * t64 - t34 * t94;
t68 = t10 * pkin(3) + t69;
t9 = -t24 * t74 + t25 * t61 + t66 * t72;
t67 = t9 * pkin(11) + t68;
t60 = sin(qJ(4));
t59 = -pkin(12) - qJ(5);
t56 = cos(pkin(13));
t50 = pkin(13) + qJ(6);
t45 = cos(t50);
t43 = sin(t50);
t42 = t56 * pkin(5) + pkin(4);
t6 = t15 * t99 + t82 * t60;
t5 = t15 * t60 - t82 * t99;
t4 = t12 * t99 + t80 * t60;
t3 = t12 * t60 - t80 * t99;
t2 = t10 * t99 + t81 * t60;
t1 = t10 * t60 - t81 * t99;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t66, -t63, 0, 0; t63, t66, 0, 0; 0, 0, 1, t91; 0, 0, 0, 1; t27, t26, t95, t87; t25, t24, -t94, t79; t36, t32, t58, t88; 0, 0, 0, 1; t12, -t11, t80, t77; t10, -t9, t81, t69; t15, -t14, t82, t78; 0, 0, 0, 1; t4, -t3, t11, t70; t2, -t1, t9, t67; t6, -t5, t14, t71; 0, 0, 0, 1; t11 * t53 + t4 * t56, t11 * t56 - t4 * t53, t3, t4 * pkin(4) + t3 * qJ(5) + t70; t2 * t56 + t9 * t53, -t2 * t53 + t9 * t56, t1, t2 * pkin(4) + t1 * qJ(5) + t67; t14 * t53 + t6 * t56, t14 * t56 - t6 * t53, t5, t6 * pkin(4) + t5 * qJ(5) + t71; 0, 0, 0, 1; t11 * t43 + t4 * t45, t11 * t45 - t4 * t43, t3, t89 * t11 - t3 * t59 + t4 * t42 + t75; t2 * t45 + t9 * t43, -t2 * t43 + t9 * t45, t1, -t1 * t59 + t2 * t42 + t89 * t9 + t68; t14 * t43 + t6 * t45, t14 * t45 - t6 * t43, t5, t89 * t14 + t6 * t42 - t5 * t59 + t76; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
