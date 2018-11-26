% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2018-11-23 15:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRRPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:26:42
% EndTime: 2018-11-23 15:26:42
% DurationCPUTime: 0.30s
% Computational Cost: add. (1589->88), mult. (1634->100), div. (0->0), fcn. (1691->22), ass. (0->74)
t51 = pkin(6) - qJ(2);
t44 = cos(t51) / 0.2e1;
t50 = pkin(6) + qJ(2);
t46 = cos(t50);
t38 = t44 + t46 / 0.2e1;
t52 = sin(pkin(12));
t55 = cos(pkin(12));
t61 = sin(qJ(2));
t29 = -t52 * t38 - t55 * t61;
t53 = sin(pkin(7));
t56 = cos(pkin(7));
t54 = sin(pkin(6));
t93 = t52 * t54;
t78 = -t29 * t53 + t56 * t93;
t43 = sin(t50) / 0.2e1;
t45 = sin(t51);
t35 = t43 + t45 / 0.2e1;
t57 = cos(pkin(6));
t80 = -t35 * t53 + t57 * t56;
t100 = pkin(5) + pkin(10);
t27 = t55 * t38 - t52 * t61;
t36 = t43 - t45 / 0.2e1;
t64 = cos(qJ(2));
t28 = t55 * t36 + t52 * t64;
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
t12 = -t27 * t72 + t28 * t60 + t55 * t70;
t99 = t12 * pkin(10);
t30 = -t52 * t36 + t55 * t64;
t14 = -t29 * t72 + t30 * t60 - t52 * t70;
t98 = t14 * pkin(10);
t39 = t44 - t46 / 0.2e1;
t17 = -t35 * t72 + t39 * t60 - t57 * t71;
t97 = t17 * pkin(10);
t96 = cos(qJ(4));
t92 = t55 * t54;
t87 = qJ(1) + 0;
t86 = t55 * pkin(1) + pkin(8) * t93 + 0;
t85 = t57 * pkin(8) + t87;
t79 = -t27 * t53 - t56 * t92;
t77 = t52 * pkin(1) - pkin(8) * t92 + 0;
t76 = t30 * pkin(2) + t78 * pkin(9) + t86;
t75 = t39 * pkin(2) + t80 * pkin(9) + t85;
t34 = t81 - t83 / 0.2e1;
t37 = t82 - t84 / 0.2e1;
t63 = cos(qJ(3));
t15 = t29 * t34 + t30 * t63 + t37 * t93;
t74 = t15 * pkin(3) + t76;
t18 = t35 * t34 + t57 * t37 + t39 * t63;
t73 = t18 * pkin(3) + t75;
t59 = sin(qJ(4));
t5 = t15 * t59 - t78 * t96;
t6 = t15 * t96 + t78 * t59;
t69 = t6 * pkin(4) + t5 * qJ(5) + t74;
t8 = t18 * t59 - t80 * t96;
t9 = t18 * t96 + t80 * t59;
t68 = t9 * pkin(4) + t8 * qJ(5) + t73;
t67 = t28 * pkin(2) + t79 * pkin(9) + t77;
t13 = t27 * t34 + t28 * t63 - t37 * t92;
t66 = t13 * pkin(3) + t67;
t3 = t13 * t59 - t79 * t96;
t4 = t13 * t96 + t79 * t59;
t65 = t4 * pkin(4) + t3 * qJ(5) + t66;
t62 = cos(qJ(6));
t58 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t55, -t52, 0, 0; t52, t55, 0, 0; 0, 0, 1, t87; 0, 0, 0, 1; t30, t29, t93, t86; t28, t27, -t92, t77; t39, t35, t57, t85; 0, 0, 0, 1; t15, -t14, t78, t76; t13, -t12, t79, t67; t18, -t17, t80, t75; 0, 0, 0, 1; t6, -t5, t14, t74 + t98; t4, -t3, t12, t66 + t99; t9, -t8, t17, t73 + t97; 0, 0, 0, 1; t14, -t6, t5, t69 + t98; t12, -t4, t3, t65 + t99; t17, -t9, t8, t68 + t97; 0, 0, 0, 1; t14 * t62 + t5 * t58, -t14 * t58 + t5 * t62, t6, t6 * pkin(11) + t100 * t14 + t69; t12 * t62 + t3 * t58, -t12 * t58 + t3 * t62, t4, t4 * pkin(11) + t100 * t12 + t65; t17 * t62 + t8 * t58, -t17 * t58 + t8 * t62, t9, t9 * pkin(11) + t100 * t17 + t68; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
