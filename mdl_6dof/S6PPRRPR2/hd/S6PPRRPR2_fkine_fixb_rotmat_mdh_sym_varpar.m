% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PPRRPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:50:16
% EndTime: 2018-11-23 14:50:16
% DurationCPUTime: 0.30s
% Computational Cost: add. (1589->88), mult. (1634->101), div. (0->0), fcn. (1691->22), ass. (0->75)
t51 = pkin(6) - pkin(12);
t44 = cos(t51) / 0.2e1;
t50 = pkin(6) + pkin(12);
t46 = cos(t50);
t36 = t44 + t46 / 0.2e1;
t52 = sin(pkin(12));
t53 = sin(pkin(11));
t57 = cos(pkin(11));
t29 = -t53 * t36 - t57 * t52;
t54 = sin(pkin(7));
t58 = cos(pkin(7));
t55 = sin(pkin(6));
t94 = t53 * t55;
t78 = -t29 * t54 + t58 * t94;
t43 = sin(t50) / 0.2e1;
t45 = sin(t51);
t34 = t43 + t45 / 0.2e1;
t59 = cos(pkin(6));
t80 = -t34 * t54 + t59 * t58;
t101 = pkin(5) + pkin(9);
t27 = t57 * t36 - t53 * t52;
t35 = t43 - t45 / 0.2e1;
t56 = cos(pkin(12));
t28 = t57 * t35 + t53 * t56;
t62 = sin(qJ(3));
t89 = pkin(7) + qJ(3);
t81 = sin(t89) / 0.2e1;
t90 = pkin(7) - qJ(3);
t83 = sin(t90);
t71 = t81 + t83 / 0.2e1;
t70 = t55 * t71;
t82 = cos(t90) / 0.2e1;
t84 = cos(t89);
t72 = t82 + t84 / 0.2e1;
t12 = -t27 * t72 + t28 * t62 + t57 * t70;
t100 = t12 * pkin(9);
t30 = -t53 * t35 + t57 * t56;
t14 = -t29 * t72 + t30 * t62 - t53 * t70;
t99 = t14 * pkin(9);
t37 = t44 - t46 / 0.2e1;
t17 = -t34 * t72 + t37 * t62 - t59 * t71;
t98 = t17 * pkin(9);
t97 = cos(qJ(4));
t93 = t57 * t55;
t91 = qJ(2) * t55;
t87 = qJ(1) + 0;
t86 = t57 * pkin(1) + t53 * t91 + 0;
t85 = t59 * qJ(2) + t87;
t79 = -t27 * t54 - t58 * t93;
t77 = t53 * pkin(1) - t57 * t91 + 0;
t76 = t30 * pkin(2) + t78 * pkin(8) + t86;
t75 = t37 * pkin(2) + t80 * pkin(8) + t85;
t38 = t81 - t83 / 0.2e1;
t39 = t82 - t84 / 0.2e1;
t64 = cos(qJ(3));
t15 = t29 * t38 + t30 * t64 + t39 * t94;
t74 = t15 * pkin(3) + t76;
t18 = t34 * t38 + t37 * t64 + t59 * t39;
t73 = t18 * pkin(3) + t75;
t61 = sin(qJ(4));
t5 = t15 * t61 - t78 * t97;
t6 = t15 * t97 + t78 * t61;
t69 = t6 * pkin(4) + t5 * qJ(5) + t74;
t8 = t18 * t61 - t80 * t97;
t9 = t18 * t97 + t80 * t61;
t68 = t9 * pkin(4) + t8 * qJ(5) + t73;
t67 = t28 * pkin(2) + t79 * pkin(8) + t77;
t13 = t27 * t38 + t28 * t64 - t39 * t93;
t66 = t13 * pkin(3) + t67;
t3 = t13 * t61 - t79 * t97;
t4 = t13 * t97 + t79 * t61;
t65 = t4 * pkin(4) + t3 * qJ(5) + t66;
t63 = cos(qJ(6));
t60 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t57, -t53, 0, 0; t53, t57, 0, 0; 0, 0, 1, t87; 0, 0, 0, 1; t30, t29, t94, t86; t28, t27, -t93, t77; t37, t34, t59, t85; 0, 0, 0, 1; t15, -t14, t78, t76; t13, -t12, t79, t67; t18, -t17, t80, t75; 0, 0, 0, 1; t6, -t5, t14, t74 + t99; t4, -t3, t12, t66 + t100; t9, -t8, t17, t73 + t98; 0, 0, 0, 1; t14, -t6, t5, t69 + t99; t12, -t4, t3, t65 + t100; t17, -t9, t8, t68 + t98; 0, 0, 0, 1; t14 * t63 + t5 * t60, -t14 * t60 + t5 * t63, t6, t6 * pkin(10) + t101 * t14 + t69; t12 * t63 + t3 * t60, -t12 * t60 + t3 * t63, t4, t4 * pkin(10) + t101 * t12 + t65; t17 * t63 + t8 * t60, -t17 * t60 + t8 * t63, t9, t9 * pkin(10) + t101 * t17 + t68; 0, 0, 0, 1;];
T_ges = t1;
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
