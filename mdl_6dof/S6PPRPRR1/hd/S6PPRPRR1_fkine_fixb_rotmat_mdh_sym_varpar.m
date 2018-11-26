% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2018-11-23 14:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:49:05
% EndTime: 2018-11-23 14:49:06
% DurationCPUTime: 0.34s
% Computational Cost: add. (1517->108), mult. (1294->133), div. (0->0), fcn. (1306->28), ass. (0->86)
t71 = pkin(7) + qJ(3);
t54 = sin(t71) / 0.2e1;
t72 = pkin(7) - qJ(3);
t61 = sin(t72);
t48 = t54 - t61 / 0.2e1;
t55 = cos(t72) / 0.2e1;
t63 = cos(t71);
t50 = t55 - t63 / 0.2e1;
t74 = sin(pkin(11));
t76 = sin(pkin(6));
t110 = t74 * t76;
t78 = cos(pkin(11));
t109 = t78 * t76;
t108 = qJ(2) * t76;
t70 = qJ(3) + pkin(13);
t107 = t74 * pkin(1) + 0;
t106 = qJ(1) + 0;
t105 = t78 * pkin(1) + t74 * t108 + 0;
t80 = cos(pkin(6));
t104 = t80 * qJ(2) + t106;
t103 = pkin(7) - t70;
t102 = pkin(7) + t70;
t68 = pkin(6) + pkin(12);
t52 = sin(t68) / 0.2e1;
t69 = pkin(6) - pkin(12);
t57 = sin(t69);
t43 = t52 + t57 / 0.2e1;
t75 = sin(pkin(7));
t79 = cos(pkin(7));
t35 = -t43 * t75 + t80 * t79;
t101 = cos(t102);
t100 = sin(t103);
t53 = cos(t69) / 0.2e1;
t58 = cos(t68);
t45 = t53 + t58 / 0.2e1;
t73 = sin(pkin(12));
t31 = t78 * t45 - t74 * t73;
t22 = -t79 * t109 - t31 * t75;
t33 = -t74 * t45 - t78 * t73;
t23 = t79 * t110 - t33 * t75;
t99 = -t78 * t108 + t107;
t98 = cos(t103) / 0.2e1;
t97 = sin(t102) / 0.2e1;
t44 = t52 - t57 / 0.2e1;
t77 = cos(pkin(12));
t34 = -t74 * t44 + t78 * t77;
t81 = pkin(8) + qJ(4);
t37 = t48 * pkin(3) - t75 * t81;
t38 = t50 * pkin(3) + t79 * t81;
t87 = cos(qJ(3));
t56 = t87 * pkin(3) + pkin(2);
t96 = t38 * t110 + t33 * t37 + t34 * t56 + t105;
t46 = t53 - t58 / 0.2e1;
t95 = t43 * t37 + t80 * t38 + t46 * t56 + t104;
t59 = sin(t70);
t89 = t100 / 0.2e1 + t97;
t88 = t76 * t89;
t90 = t101 / 0.2e1 + t98;
t11 = -t33 * t90 + t34 * t59 - t74 * t88;
t41 = t97 - t100 / 0.2e1;
t42 = t98 - t101 / 0.2e1;
t62 = cos(t70);
t12 = t42 * t110 + t33 * t41 + t34 * t62;
t94 = t12 * pkin(4) + t11 * pkin(9) + t96;
t32 = t78 * t44 + t74 * t77;
t93 = t31 * t37 + t32 * t56 + (-qJ(2) - t38) * t109 + t107;
t14 = -t43 * t90 + t46 * t59 - t80 * t89;
t15 = t43 * t41 + t80 * t42 + t46 * t62;
t92 = t15 * pkin(4) + t14 * pkin(9) + t95;
t10 = -t42 * t109 + t31 * t41 + t32 * t62;
t9 = -t31 * t90 + t32 * t59 + t78 * t88;
t91 = t10 * pkin(4) + t9 * pkin(9) + t93;
t86 = cos(qJ(5));
t85 = cos(qJ(6));
t84 = sin(qJ(3));
t83 = sin(qJ(5));
t82 = sin(qJ(6));
t49 = t55 + t63 / 0.2e1;
t47 = t54 + t61 / 0.2e1;
t6 = t15 * t86 + t35 * t83;
t5 = t15 * t83 - t35 * t86;
t4 = t12 * t86 + t23 * t83;
t3 = t12 * t83 - t23 * t86;
t2 = t10 * t86 + t22 * t83;
t1 = t10 * t83 - t22 * t86;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t78, -t74, 0, 0; t74, t78, 0, 0; 0, 0, 1, t106; 0, 0, 0, 1; t34, t33, t110, t105; t32, t31, -t109, t99; t46, t43, t80, t104; 0, 0, 0, 1; t50 * t110 + t33 * t48 + t34 * t87, t47 * t110 + t33 * t49 - t34 * t84, t23, t34 * pkin(2) + t23 * pkin(8) + t105; -t50 * t109 + t31 * t48 + t32 * t87, -t47 * t109 + t31 * t49 - t32 * t84, t22, t32 * pkin(2) + t22 * pkin(8) + t99; t43 * t48 + t46 * t87 + t80 * t50, t43 * t49 - t46 * t84 + t80 * t47, t35, t46 * pkin(2) + t35 * pkin(8) + t104; 0, 0, 0, 1; t12, -t11, t23, t96; t10, -t9, t22, t93; t15, -t14, t35, t95; 0, 0, 0, 1; t4, -t3, t11, t94; t2, -t1, t9, t91; t6, -t5, t14, t92; 0, 0, 0, 1; t11 * t82 + t4 * t85, t11 * t85 - t4 * t82, t3, t4 * pkin(5) + t3 * pkin(10) + t94; t2 * t85 + t9 * t82, -t2 * t82 + t9 * t85, t1, t2 * pkin(5) + t1 * pkin(10) + t91; t14 * t82 + t6 * t85, t14 * t85 - t6 * t82, t5, t6 * pkin(5) + t5 * pkin(10) + t92; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
