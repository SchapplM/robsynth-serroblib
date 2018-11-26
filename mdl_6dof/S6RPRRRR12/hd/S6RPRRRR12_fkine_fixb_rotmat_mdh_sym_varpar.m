% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRRR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:40:40
% EndTime: 2018-11-23 16:40:40
% DurationCPUTime: 0.44s
% Computational Cost: add. (3259->102), mult. (3274->130), div. (0->0), fcn. (3324->30), ass. (0->93)
t81 = sin(pkin(6));
t90 = sin(qJ(1));
t120 = t90 * t81;
t75 = pkin(6) - pkin(14);
t64 = cos(t75) / 0.2e1;
t74 = pkin(6) + pkin(14);
t68 = cos(t74);
t52 = t64 + t68 / 0.2e1;
t78 = sin(pkin(14));
t95 = cos(qJ(1));
t47 = -t90 * t52 - t95 * t78;
t80 = sin(pkin(7));
t84 = cos(pkin(7));
t39 = t84 * t120 - t47 * t80;
t63 = sin(t74) / 0.2e1;
t67 = sin(t75);
t50 = t63 + t67 / 0.2e1;
t85 = cos(pkin(6));
t44 = -t50 * t80 + t85 * t84;
t53 = t64 - t68 / 0.2e1;
t76 = pkin(7) + qJ(3);
t65 = sin(t76) / 0.2e1;
t77 = pkin(7) - qJ(3);
t69 = sin(t77);
t55 = t65 + t69 / 0.2e1;
t66 = cos(t77) / 0.2e1;
t70 = cos(t76);
t58 = t66 + t70 / 0.2e1;
t89 = sin(qJ(3));
t32 = t50 * t58 - t53 * t89 + t85 * t55;
t79 = sin(pkin(8));
t83 = cos(pkin(8));
t23 = -t32 * t79 + t44 * t83;
t51 = t63 - t67 / 0.2e1;
t82 = cos(pkin(14));
t48 = -t90 * t51 + t95 * t82;
t28 = t120 * t55 + t47 * t58 - t48 * t89;
t19 = -t28 * t79 + t39 * t83;
t119 = t95 * t81;
t45 = t95 * t52 - t90 * t78;
t46 = t95 * t51 + t90 * t82;
t26 = -t119 * t55 + t45 * t58 - t46 * t89;
t38 = -t119 * t84 - t45 * t80;
t18 = -t26 * t79 + t38 * t83;
t118 = qJ(2) * t81;
t117 = pkin(8) - qJ(4);
t116 = pkin(8) + qJ(4);
t115 = pkin(9) + 0;
t113 = t85 * qJ(2) + t115;
t112 = t95 * pkin(1) + t90 * t118 + 0;
t111 = cos(t116);
t110 = sin(t117);
t109 = cos(t117) / 0.2e1;
t108 = sin(t116) / 0.2e1;
t107 = t90 * pkin(1) - t118 * t95 + 0;
t106 = t53 * pkin(2) + t44 * pkin(10) + t113;
t105 = t48 * pkin(2) + t39 * pkin(10) + t112;
t104 = t109 + t111 / 0.2e1;
t103 = t108 + t110 / 0.2e1;
t56 = t65 - t69 / 0.2e1;
t59 = t66 - t70 / 0.2e1;
t94 = cos(qJ(3));
t33 = t50 * t56 + t53 * t94 + t85 * t59;
t102 = t33 * pkin(3) + t23 * pkin(11) + t106;
t29 = t120 * t59 + t47 * t56 + t48 * t94;
t101 = t29 * pkin(3) + t19 * pkin(11) + t105;
t100 = t46 * pkin(2) + pkin(10) * t38 + t107;
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
t27 = -t119 * t59 + t45 * t56 + t46 * t94;
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
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t95, -t90, 0, 0; t90, t95, 0, 0; 0, 0, 1, t115; 0, 0, 0, 1; t48, t47, t120, t112; t46, t45, -t119, t107; t53, t50, t85, t113; 0, 0, 0, 1; t29, t28, t39, t105; t27, t26, t38, t100; t33, t32, t44, t106; 0, 0, 0, 1; t12, -t11, t19, t101; t10, -t9, t18, t97; t15, -t14, t23, t102; 0, 0, 0, 1; t4, -t3, t11, t98; t2, -t1, t9, t96; t6, -t5, t14, t99; 0, 0, 0, 1; t11 * t86 + t4 * t91, t11 * t91 - t4 * t86, t3, t4 * pkin(5) + t3 * pkin(13) + t98; t2 * t91 + t9 * t86, -t2 * t86 + t9 * t91, t1, t2 * pkin(5) + t1 * pkin(13) + t96; t14 * t86 + t6 * t91, t14 * t91 - t6 * t86, t5, t6 * pkin(5) + t5 * pkin(13) + t99; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
