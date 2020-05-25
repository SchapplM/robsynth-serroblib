% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2018-11-23 16:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:07:59
% EndTime: 2018-11-23 16:07:59
% DurationCPUTime: 0.30s
% Computational Cost: add. (1517->108), mult. (1294->133), div. (0->0), fcn. (1306->28), ass. (0->86)
t71 = pkin(7) + qJ(3);
t54 = sin(t71) / 0.2e1;
t72 = pkin(7) - qJ(3);
t61 = sin(t72);
t48 = t54 - t61 / 0.2e1;
t55 = cos(t72) / 0.2e1;
t63 = cos(t71);
t50 = t55 - t63 / 0.2e1;
t75 = sin(pkin(6));
t83 = sin(qJ(1));
t110 = t83 * t75;
t87 = cos(qJ(1));
t109 = t87 * t75;
t108 = qJ(2) * t75;
t107 = pkin(8) + 0;
t70 = qJ(3) + pkin(13);
t106 = t83 * pkin(1) + 0;
t78 = cos(pkin(6));
t105 = t78 * qJ(2) + t107;
t104 = t87 * pkin(1) + t83 * t108 + 0;
t103 = pkin(7) - t70;
t102 = pkin(7) + t70;
t68 = pkin(6) + pkin(12);
t52 = sin(t68) / 0.2e1;
t69 = pkin(6) - pkin(12);
t57 = sin(t69);
t43 = t52 + t57 / 0.2e1;
t74 = sin(pkin(7));
t77 = cos(pkin(7));
t31 = -t43 * t74 + t78 * t77;
t101 = cos(t102);
t100 = sin(t103);
t53 = cos(t69) / 0.2e1;
t58 = cos(t68);
t45 = t53 + t58 / 0.2e1;
t73 = sin(pkin(12));
t32 = t87 * t45 - t83 * t73;
t22 = -t77 * t109 - t32 * t74;
t34 = -t83 * t45 - t87 * t73;
t23 = t77 * t110 - t34 * t74;
t79 = pkin(9) + qJ(4);
t36 = t48 * pkin(3) - t74 * t79;
t37 = t50 * pkin(3) + t77 * t79;
t46 = t53 - t58 / 0.2e1;
t86 = cos(qJ(3));
t56 = t86 * pkin(3) + pkin(2);
t99 = t43 * t36 + t78 * t37 + t46 * t56 + t105;
t98 = -t87 * t108 + t106;
t97 = cos(t103) / 0.2e1;
t96 = sin(t102) / 0.2e1;
t44 = t52 - t57 / 0.2e1;
t76 = cos(pkin(12));
t35 = -t83 * t44 + t87 * t76;
t95 = t37 * t110 + t34 * t36 + t35 * t56 + t104;
t59 = sin(t70);
t89 = t100 / 0.2e1 + t96;
t90 = t101 / 0.2e1 + t97;
t14 = -t43 * t90 + t46 * t59 - t78 * t89;
t41 = t96 - t100 / 0.2e1;
t42 = t97 - t101 / 0.2e1;
t62 = cos(t70);
t15 = t43 * t41 + t78 * t42 + t46 * t62;
t94 = t15 * pkin(4) + t14 * pkin(10) + t99;
t88 = t75 * t89;
t11 = -t34 * t90 + t35 * t59 - t83 * t88;
t12 = t42 * t110 + t34 * t41 + t35 * t62;
t93 = t12 * pkin(4) + t11 * pkin(10) + t95;
t33 = t87 * t44 + t83 * t76;
t92 = t32 * t36 + t33 * t56 + (-qJ(2) - t37) * t109 + t106;
t10 = -t42 * t109 + t32 * t41 + t33 * t62;
t9 = -t32 * t90 + t33 * t59 + t87 * t88;
t91 = t10 * pkin(4) + t9 * pkin(10) + t92;
t85 = cos(qJ(5));
t84 = cos(qJ(6));
t82 = sin(qJ(3));
t81 = sin(qJ(5));
t80 = sin(qJ(6));
t49 = t55 + t63 / 0.2e1;
t47 = t54 + t61 / 0.2e1;
t6 = t15 * t85 + t31 * t81;
t5 = t15 * t81 - t31 * t85;
t4 = t12 * t85 + t23 * t81;
t3 = t12 * t81 - t23 * t85;
t2 = t10 * t85 + t22 * t81;
t1 = t10 * t81 - t22 * t85;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t87, -t83, 0, 0; t83, t87, 0, 0; 0, 0, 1, t107; 0, 0, 0, 1; t35, t34, t110, t104; t33, t32, -t109, t98; t46, t43, t78, t105; 0, 0, 0, 1; t50 * t110 + t34 * t48 + t35 * t86, t47 * t110 + t34 * t49 - t35 * t82, t23, t35 * pkin(2) + t23 * pkin(9) + t104; -t50 * t109 + t32 * t48 + t33 * t86, -t47 * t109 + t32 * t49 - t33 * t82, t22, t33 * pkin(2) + t22 * pkin(9) + t98; t43 * t48 + t46 * t86 + t78 * t50, t43 * t49 - t46 * t82 + t78 * t47, t31, t46 * pkin(2) + t31 * pkin(9) + t105; 0, 0, 0, 1; t12, -t11, t23, t95; t10, -t9, t22, t92; t15, -t14, t31, t99; 0, 0, 0, 1; t4, -t3, t11, t93; t2, -t1, t9, t91; t6, -t5, t14, t94; 0, 0, 0, 1; t11 * t80 + t4 * t84, t11 * t84 - t4 * t80, t3, t4 * pkin(5) + t3 * pkin(11) + t93; t2 * t84 + t9 * t80, -t2 * t80 + t9 * t84, t1, t2 * pkin(5) + t1 * pkin(11) + t91; t14 * t80 + t6 * t84, t14 * t84 - t6 * t80, t5, t6 * pkin(5) + t5 * pkin(11) + t94; 0, 0, 0, 1;];
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
