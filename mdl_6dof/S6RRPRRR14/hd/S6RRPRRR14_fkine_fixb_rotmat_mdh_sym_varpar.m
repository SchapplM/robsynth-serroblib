% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRR14_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:09:41
% EndTime: 2018-12-10 18:09:41
% DurationCPUTime: 0.46s
% Computational Cost: add. (3259->104), mult. (3274->132), div. (0->0), fcn. (3324->30), ass. (0->95)
t74 = pkin(7) + pkin(14);
t63 = sin(t74) / 0.2e1;
t75 = pkin(7) - pkin(14);
t67 = sin(t75);
t50 = t63 + t67 / 0.2e1;
t64 = cos(t75) / 0.2e1;
t68 = cos(t74);
t52 = t64 + t68 / 0.2e1;
t76 = pkin(6) + qJ(2);
t65 = sin(t76) / 0.2e1;
t77 = pkin(6) - qJ(2);
t69 = sin(t77);
t55 = t65 + t69 / 0.2e1;
t66 = cos(t77) / 0.2e1;
t70 = cos(t76);
t59 = t66 - t70 / 0.2e1;
t78 = sin(pkin(14));
t85 = cos(pkin(6));
t32 = t85 * t50 + t55 * t52 - t59 * t78;
t84 = cos(pkin(7));
t121 = t85 * t84;
t80 = sin(pkin(7));
t44 = -t55 * t80 + t121;
t79 = sin(pkin(8));
t83 = cos(pkin(8));
t23 = -t32 * t79 + t44 * t83;
t81 = sin(pkin(6));
t90 = sin(qJ(1));
t120 = t90 * t81;
t58 = t66 + t70 / 0.2e1;
t89 = sin(qJ(2));
t95 = cos(qJ(1));
t47 = -t90 * t58 - t95 * t89;
t56 = t65 - t69 / 0.2e1;
t94 = cos(qJ(2));
t48 = -t90 * t56 + t95 * t94;
t28 = t120 * t50 + t47 * t52 - t48 * t78;
t114 = t84 * t120;
t39 = -t47 * t80 + t114;
t19 = -t28 * t79 + t39 * t83;
t119 = t95 * t81;
t45 = t95 * t58 - t90 * t89;
t46 = t95 * t56 + t90 * t94;
t26 = -t119 * t50 + t45 * t52 - t46 * t78;
t38 = -t119 * t84 - t45 * t80;
t18 = -t26 * t79 + t38 * t83;
t118 = qJ(3) * t80;
t117 = pkin(8) - qJ(4);
t116 = pkin(8) + qJ(4);
t115 = pkin(9) + 0;
t113 = t85 * pkin(10) + t115;
t112 = t95 * pkin(1) + pkin(10) * t120 + 0;
t111 = cos(t116);
t110 = sin(t117);
t109 = cos(t117) / 0.2e1;
t108 = sin(t116) / 0.2e1;
t107 = t90 * pkin(1) - pkin(10) * t119 + 0;
t106 = t59 * pkin(2) + qJ(3) * t121 - t118 * t55 + t113;
t105 = t48 * pkin(2) + qJ(3) * t114 - t118 * t47 + t112;
t104 = t109 + t111 / 0.2e1;
t103 = t108 + t110 / 0.2e1;
t51 = t63 - t67 / 0.2e1;
t53 = t64 - t68 / 0.2e1;
t82 = cos(pkin(14));
t33 = t55 * t51 + t85 * t53 + t59 * t82;
t102 = t33 * pkin(3) + t23 * pkin(11) + t106;
t29 = t120 * t53 + t47 * t51 + t48 * t82;
t101 = t29 * pkin(3) + t19 * pkin(11) + t105;
t100 = t46 * pkin(2) + qJ(3) * t38 + t107;
t88 = sin(qJ(4));
t14 = -t103 * t44 - t104 * t32 + t33 * t88;
t54 = t108 - t110 / 0.2e1;
t57 = t109 - t111 / 0.2e1;
t93 = cos(qJ(4));
t15 = t32 * t54 + t33 * t93 + t44 * t57;
t99 = t15 * pkin(4) + t14 * pkin(12) + t102;
t11 = -t103 * t39 - t104 * t28 + t29 * t88;
t12 = t28 * t54 + t29 * t93 + t39 * t57;
t98 = t12 * pkin(4) + t11 * pkin(12) + t101;
t27 = -t119 * t53 + t45 * t51 + t46 * t82;
t97 = t27 * pkin(3) + t18 * pkin(11) + t100;
t10 = t26 * t54 + t27 * t93 + t38 * t57;
t9 = -t103 * t38 - t104 * t26 + t27 * t88;
t96 = t10 * pkin(4) + t9 * pkin(12) + t97;
t92 = cos(qJ(5));
t91 = cos(qJ(6));
t87 = sin(qJ(5));
t86 = sin(qJ(6));
t6 = t15 * t92 + t23 * t87;
t5 = t15 * t87 - t23 * t92;
t4 = t12 * t92 + t19 * t87;
t3 = t12 * t87 - t19 * t92;
t2 = t10 * t92 + t18 * t87;
t1 = t10 * t87 - t18 * t92;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t95, -t90, 0, 0; t90, t95, 0, 0; 0, 0, 1, t115; 0, 0, 0, 1; t48, t47, t120, t112; t46, t45, -t119, t107; t59, t55, t85, t113; 0, 0, 0, 1; t29, t28, t39, t105; t27, t26, t38, t100; t33, t32, t44, t106; 0, 0, 0, 1; t12, -t11, t19, t101; t10, -t9, t18, t97; t15, -t14, t23, t102; 0, 0, 0, 1; t4, -t3, t11, t98; t2, -t1, t9, t96; t6, -t5, t14, t99; 0, 0, 0, 1; t11 * t86 + t4 * t91, t11 * t91 - t4 * t86, t3, t4 * pkin(5) + t3 * pkin(13) + t98; t2 * t91 + t9 * t86, -t2 * t86 + t9 * t91, t1, t2 * pkin(5) + t1 * pkin(13) + t96; t14 * t86 + t6 * t91, t14 * t91 - t6 * t86, t5, t6 * pkin(5) + t5 * pkin(13) + t99; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
