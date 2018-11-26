% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2018-11-23 16:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:23:03
% EndTime: 2018-11-23 16:23:03
% DurationCPUTime: 0.28s
% Computational Cost: add. (1589->88), mult. (1634->101), div. (0->0), fcn. (1691->22), ass. (0->75)
t51 = pkin(6) - pkin(12);
t44 = cos(t51) / 0.2e1;
t50 = pkin(6) + pkin(12);
t46 = cos(t50);
t36 = t44 + t46 / 0.2e1;
t52 = sin(pkin(12));
t61 = sin(qJ(1));
t64 = cos(qJ(1));
t29 = -t61 * t36 - t64 * t52;
t53 = sin(pkin(7));
t56 = cos(pkin(7));
t54 = sin(pkin(6));
t93 = t61 * t54;
t78 = -t29 * t53 + t56 * t93;
t43 = sin(t50) / 0.2e1;
t45 = sin(t51);
t34 = t43 + t45 / 0.2e1;
t57 = cos(pkin(6));
t80 = -t34 * t53 + t57 * t56;
t101 = pkin(5) + pkin(10);
t27 = t64 * t36 - t61 * t52;
t35 = t43 - t45 / 0.2e1;
t55 = cos(pkin(12));
t28 = t64 * t35 + t61 * t55;
t60 = sin(qJ(3));
t89 = pkin(7) + qJ(3);
t81 = sin(t89) / 0.2e1;
t90 = pkin(7) - qJ(3);
t83 = sin(t90);
t71 = t81 + t83 / 0.2e1;
t70 = t54 * t71;
t82 = cos(t90) / 0.2e1;
t84 = cos(t89);
t72 = t82 + t84 / 0.2e1;
t12 = -t27 * t72 + t28 * t60 + t64 * t70;
t100 = t12 * pkin(10);
t30 = -t61 * t35 + t64 * t55;
t14 = -t29 * t72 + t30 * t60 - t61 * t70;
t99 = t14 * pkin(10);
t37 = t44 - t46 / 0.2e1;
t17 = -t34 * t72 + t37 * t60 - t57 * t71;
t98 = t17 * pkin(10);
t97 = cos(qJ(4));
t92 = t64 * t54;
t91 = qJ(2) * t54;
t88 = pkin(8) + 0;
t86 = t57 * qJ(2) + t88;
t85 = t64 * pkin(1) + t61 * t91 + 0;
t79 = -t27 * t53 - t56 * t92;
t77 = t61 * pkin(1) - t64 * t91 + 0;
t76 = t37 * pkin(2) + t80 * pkin(9) + t86;
t75 = t30 * pkin(2) + t78 * pkin(9) + t85;
t38 = t81 - t83 / 0.2e1;
t39 = t82 - t84 / 0.2e1;
t63 = cos(qJ(3));
t18 = t34 * t38 + t37 * t63 + t57 * t39;
t74 = t18 * pkin(3) + t76;
t15 = t29 * t38 + t30 * t63 + t39 * t93;
t73 = t15 * pkin(3) + t75;
t59 = sin(qJ(4));
t8 = t18 * t59 - t80 * t97;
t9 = t18 * t97 + t80 * t59;
t69 = t9 * pkin(4) + t8 * qJ(5) + t74;
t5 = t15 * t59 - t78 * t97;
t6 = t15 * t97 + t78 * t59;
t68 = t6 * pkin(4) + t5 * qJ(5) + t73;
t67 = t28 * pkin(2) + t79 * pkin(9) + t77;
t13 = t27 * t38 + t28 * t63 - t39 * t92;
t66 = t13 * pkin(3) + t67;
t3 = t13 * t59 - t79 * t97;
t4 = t13 * t97 + t79 * t59;
t65 = t4 * pkin(4) + t3 * qJ(5) + t66;
t62 = cos(qJ(6));
t58 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t64, -t61, 0, 0; t61, t64, 0, 0; 0, 0, 1, t88; 0, 0, 0, 1; t30, t29, t93, t85; t28, t27, -t92, t77; t37, t34, t57, t86; 0, 0, 0, 1; t15, -t14, t78, t75; t13, -t12, t79, t67; t18, -t17, t80, t76; 0, 0, 0, 1; t6, -t5, t14, t73 + t99; t4, -t3, t12, t66 + t100; t9, -t8, t17, t74 + t98; 0, 0, 0, 1; t14, -t6, t5, t68 + t99; t12, -t4, t3, t65 + t100; t17, -t9, t8, t69 + t98; 0, 0, 0, 1; t14 * t62 + t5 * t58, -t14 * t58 + t5 * t62, t6, t6 * pkin(11) + t101 * t14 + t68; t12 * t62 + t3 * t58, -t12 * t58 + t3 * t62, t4, t4 * pkin(11) + t101 * t12 + t65; t17 * t62 + t8 * t58, -t17 * t58 + t8 * t62, t9, t9 * pkin(11) + t101 * t17 + t69; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
