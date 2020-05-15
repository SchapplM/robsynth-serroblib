% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRRPR15_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:25:00
% EndTime: 2018-11-23 18:25:00
% DurationCPUTime: 0.28s
% Computational Cost: add. (1589->88), mult. (1634->100), div. (0->0), fcn. (1691->22), ass. (0->74)
t51 = pkin(6) - qJ(2);
t44 = cos(t51) / 0.2e1;
t50 = pkin(6) + qJ(2);
t46 = cos(t50);
t38 = t44 + t46 / 0.2e1;
t59 = sin(qJ(2));
t60 = sin(qJ(1));
t64 = cos(qJ(1));
t29 = -t60 * t38 - t64 * t59;
t52 = sin(pkin(7));
t54 = cos(pkin(7));
t53 = sin(pkin(6));
t92 = t60 * t53;
t78 = -t29 * t52 + t54 * t92;
t43 = sin(t50) / 0.2e1;
t45 = sin(t51);
t35 = t43 + t45 / 0.2e1;
t55 = cos(pkin(6));
t80 = -t35 * t52 + t55 * t54;
t100 = pkin(5) + pkin(11);
t27 = t64 * t38 - t60 * t59;
t36 = t43 - t45 / 0.2e1;
t63 = cos(qJ(2));
t28 = t64 * t36 + t60 * t63;
t58 = sin(qJ(3));
t89 = pkin(7) + qJ(3);
t81 = sin(t89) / 0.2e1;
t90 = pkin(7) - qJ(3);
t83 = sin(t90);
t71 = t81 + t83 / 0.2e1;
t70 = t53 * t71;
t82 = cos(t90) / 0.2e1;
t84 = cos(t89);
t72 = t82 + t84 / 0.2e1;
t12 = -t27 * t72 + t28 * t58 + t64 * t70;
t99 = t12 * pkin(11);
t30 = -t60 * t36 + t64 * t63;
t14 = -t29 * t72 + t30 * t58 - t60 * t70;
t98 = t14 * pkin(11);
t39 = t44 - t46 / 0.2e1;
t17 = -t35 * t72 + t39 * t58 - t55 * t71;
t97 = t17 * pkin(11);
t96 = cos(qJ(4));
t91 = t64 * t53;
t88 = pkin(8) + 0;
t86 = t55 * pkin(9) + t88;
t85 = t64 * pkin(1) + pkin(9) * t92 + 0;
t79 = -t27 * t52 - t54 * t91;
t77 = t60 * pkin(1) - pkin(9) * t91 + 0;
t76 = t39 * pkin(2) + t80 * pkin(10) + t86;
t75 = t30 * pkin(2) + t78 * pkin(10) + t85;
t34 = t81 - t83 / 0.2e1;
t37 = t82 - t84 / 0.2e1;
t62 = cos(qJ(3));
t18 = t35 * t34 + t55 * t37 + t39 * t62;
t74 = t18 * pkin(3) + t76;
t15 = t29 * t34 + t30 * t62 + t37 * t92;
t73 = t15 * pkin(3) + t75;
t57 = sin(qJ(4));
t8 = t18 * t57 - t80 * t96;
t9 = t18 * t96 + t80 * t57;
t69 = t9 * pkin(4) + t8 * qJ(5) + t74;
t5 = t15 * t57 - t78 * t96;
t6 = t15 * t96 + t78 * t57;
t68 = t6 * pkin(4) + t5 * qJ(5) + t73;
t67 = t28 * pkin(2) + t79 * pkin(10) + t77;
t13 = t27 * t34 + t28 * t62 - t37 * t91;
t66 = t13 * pkin(3) + t67;
t3 = t13 * t57 - t79 * t96;
t4 = t13 * t96 + t79 * t57;
t65 = t4 * pkin(4) + t3 * qJ(5) + t66;
t61 = cos(qJ(6));
t56 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t64, -t60, 0, 0; t60, t64, 0, 0; 0, 0, 1, t88; 0, 0, 0, 1; t30, t29, t92, t85; t28, t27, -t91, t77; t39, t35, t55, t86; 0, 0, 0, 1; t15, -t14, t78, t75; t13, -t12, t79, t67; t18, -t17, t80, t76; 0, 0, 0, 1; t6, -t5, t14, t73 + t98; t4, -t3, t12, t66 + t99; t9, -t8, t17, t74 + t97; 0, 0, 0, 1; t14, -t6, t5, t68 + t98; t12, -t4, t3, t65 + t99; t17, -t9, t8, t69 + t97; 0, 0, 0, 1; t14 * t61 + t5 * t56, -t14 * t56 + t5 * t61, t6, t6 * pkin(12) + t100 * t14 + t68; t12 * t61 + t3 * t56, -t12 * t56 + t3 * t61, t4, t4 * pkin(12) + t100 * t12 + t65; t17 * t61 + t8 * t56, -t17 * t56 + t8 * t61, t9, t9 * pkin(12) + t100 * t17 + t69; 0, 0, 0, 1;];
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
