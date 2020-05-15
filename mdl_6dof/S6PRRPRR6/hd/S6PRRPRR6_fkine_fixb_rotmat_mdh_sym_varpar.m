% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRRPRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:17:48
% EndTime: 2018-11-23 15:17:49
% DurationCPUTime: 0.32s
% Computational Cost: add. (1430->93), mult. (1432->116), div. (0->0), fcn. (1486->24), ass. (0->75)
t61 = pkin(6) - qJ(2);
t50 = cos(t61) / 0.2e1;
t60 = pkin(6) + qJ(2);
t55 = cos(t60);
t44 = t50 + t55 / 0.2e1;
t63 = sin(pkin(12));
t67 = cos(pkin(12));
t73 = sin(qJ(2));
t35 = -t63 * t44 - t67 * t73;
t64 = sin(pkin(7));
t68 = cos(pkin(7));
t65 = sin(pkin(6));
t99 = t63 * t65;
t24 = -t35 * t64 + t68 * t99;
t49 = sin(t60) / 0.2e1;
t53 = sin(t61);
t41 = t49 + t53 / 0.2e1;
t69 = cos(pkin(6));
t32 = -t41 * t64 + t69 * t68;
t33 = t67 * t44 - t63 * t73;
t98 = t67 * t65;
t23 = -t33 * t64 - t68 * t98;
t62 = sin(pkin(13));
t104 = t23 * t62;
t103 = t24 * t62;
t102 = t32 * t62;
t96 = pkin(7) - qJ(3);
t95 = pkin(7) + qJ(3);
t93 = qJ(1) + 0;
t92 = t67 * pkin(1) + pkin(8) * t99 + 0;
t91 = t69 * pkin(8) + t93;
t90 = cos(t95);
t89 = sin(t96);
t88 = cos(t96) / 0.2e1;
t87 = sin(t95) / 0.2e1;
t86 = t63 * pkin(1) - pkin(8) * t98 + 0;
t42 = t49 - t53 / 0.2e1;
t76 = cos(qJ(2));
t36 = -t63 * t42 + t67 * t76;
t85 = t36 * pkin(2) + t24 * pkin(9) + t92;
t45 = t50 - t55 / 0.2e1;
t84 = t45 * pkin(2) + t32 * pkin(9) + t91;
t83 = t88 + t90 / 0.2e1;
t82 = t87 + t89 / 0.2e1;
t72 = sin(qJ(3));
t80 = t65 * t82;
t13 = -t35 * t83 + t36 * t72 - t63 * t80;
t40 = t87 - t89 / 0.2e1;
t43 = t88 - t90 / 0.2e1;
t75 = cos(qJ(3));
t14 = t35 * t40 + t36 * t75 + t43 * t99;
t66 = cos(pkin(13));
t51 = t66 * pkin(4) + pkin(3);
t70 = -pkin(10) - qJ(4);
t81 = pkin(4) * t103 - t13 * t70 + t14 * t51 + t85;
t17 = -t41 * t83 + t45 * t72 - t69 * t82;
t18 = t41 * t40 + t69 * t43 + t45 * t75;
t79 = pkin(4) * t102 - t17 * t70 + t18 * t51 + t84;
t34 = t67 * t42 + t63 * t76;
t78 = t34 * pkin(2) + t23 * pkin(9) + t86;
t11 = -t33 * t83 + t34 * t72 + t67 * t80;
t12 = t33 * t40 + t34 * t75 - t43 * t98;
t77 = pkin(4) * t104 - t11 * t70 + t12 * t51 + t78;
t74 = cos(qJ(6));
t71 = sin(qJ(6));
t59 = pkin(13) + qJ(5);
t54 = cos(t59);
t52 = sin(t59);
t6 = t18 * t54 + t32 * t52;
t5 = t18 * t52 - t32 * t54;
t4 = t14 * t54 + t24 * t52;
t3 = t14 * t52 - t24 * t54;
t2 = t12 * t54 + t23 * t52;
t1 = t12 * t52 - t23 * t54;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t67, -t63, 0, 0; t63, t67, 0, 0; 0, 0, 1, t93; 0, 0, 0, 1; t36, t35, t99, t92; t34, t33, -t98, t86; t45, t41, t69, t91; 0, 0, 0, 1; t14, -t13, t24, t85; t12, -t11, t23, t78; t18, -t17, t32, t84; 0, 0, 0, 1; t14 * t66 + t103, -t14 * t62 + t24 * t66, t13, t14 * pkin(3) + t13 * qJ(4) + t85; t12 * t66 + t104, -t12 * t62 + t23 * t66, t11, t12 * pkin(3) + t11 * qJ(4) + t78; t18 * t66 + t102, -t18 * t62 + t32 * t66, t17, t18 * pkin(3) + t17 * qJ(4) + t84; 0, 0, 0, 1; t4, -t3, t13, t81; t2, -t1, t11, t77; t6, -t5, t17, t79; 0, 0, 0, 1; t13 * t71 + t4 * t74, t13 * t74 - t4 * t71, t3, t4 * pkin(5) + t3 * pkin(11) + t81; t11 * t71 + t2 * t74, t11 * t74 - t2 * t71, t1, t2 * pkin(5) + t1 * pkin(11) + t77; t17 * t71 + t6 * t74, t17 * t74 - t6 * t71, t5, t6 * pkin(5) + t5 * pkin(11) + t79; 0, 0, 0, 1;];
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
