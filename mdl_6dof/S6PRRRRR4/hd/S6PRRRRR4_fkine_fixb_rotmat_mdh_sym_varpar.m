% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRRRRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:34:13
% EndTime: 2018-11-23 15:34:13
% DurationCPUTime: 0.30s
% Computational Cost: add. (1430->93), mult. (1432->116), div. (0->0), fcn. (1486->24), ass. (0->75)
t60 = pkin(6) - qJ(2);
t50 = cos(t60) / 0.2e1;
t59 = pkin(6) + qJ(2);
t53 = cos(t59);
t44 = t50 + t53 / 0.2e1;
t62 = sin(pkin(13));
t65 = cos(pkin(13));
t71 = sin(qJ(2));
t35 = -t62 * t44 - t65 * t71;
t63 = sin(pkin(7));
t66 = cos(pkin(7));
t64 = sin(pkin(6));
t99 = t62 * t64;
t24 = -t35 * t63 + t66 * t99;
t49 = sin(t59) / 0.2e1;
t52 = sin(t60);
t41 = t49 + t52 / 0.2e1;
t67 = cos(pkin(6));
t32 = -t41 * t63 + t67 * t66;
t33 = t65 * t44 - t62 * t71;
t98 = t65 * t64;
t23 = -t33 * t63 - t66 * t98;
t69 = sin(qJ(4));
t104 = t23 * t69;
t103 = t24 * t69;
t102 = t32 * t69;
t96 = pkin(7) - qJ(3);
t95 = pkin(7) + qJ(3);
t93 = qJ(1) + 0;
t92 = t65 * pkin(1) + pkin(8) * t99 + 0;
t91 = t67 * pkin(8) + t93;
t90 = cos(t95);
t89 = sin(t96);
t88 = cos(t96) / 0.2e1;
t87 = sin(t95) / 0.2e1;
t86 = t62 * pkin(1) - pkin(8) * t98 + 0;
t42 = t49 - t52 / 0.2e1;
t75 = cos(qJ(2));
t36 = -t62 * t42 + t65 * t75;
t85 = t36 * pkin(2) + t24 * pkin(9) + t92;
t45 = t50 - t53 / 0.2e1;
t84 = t45 * pkin(2) + t32 * pkin(9) + t91;
t83 = t88 + t90 / 0.2e1;
t82 = t87 + t89 / 0.2e1;
t70 = sin(qJ(3));
t80 = t64 * t82;
t13 = -t35 * t83 + t36 * t70 - t62 * t80;
t40 = t87 - t89 / 0.2e1;
t43 = t88 - t90 / 0.2e1;
t74 = cos(qJ(3));
t14 = t35 * t40 + t36 * t74 + t43 * t99;
t73 = cos(qJ(4));
t51 = t73 * pkin(4) + pkin(3);
t76 = -pkin(11) - pkin(10);
t81 = pkin(4) * t103 - t13 * t76 + t14 * t51 + t85;
t17 = -t41 * t83 + t45 * t70 - t67 * t82;
t18 = t41 * t40 + t67 * t43 + t45 * t74;
t79 = pkin(4) * t102 - t17 * t76 + t18 * t51 + t84;
t34 = t65 * t42 + t62 * t75;
t78 = t34 * pkin(2) + t23 * pkin(9) + t86;
t11 = -t33 * t83 + t34 * t70 + t65 * t80;
t12 = t33 * t40 + t34 * t74 - t43 * t98;
t77 = pkin(4) * t104 - t11 * t76 + t12 * t51 + t78;
t72 = cos(qJ(6));
t68 = sin(qJ(6));
t61 = qJ(4) + qJ(5);
t55 = cos(t61);
t54 = sin(t61);
t6 = t18 * t55 + t32 * t54;
t5 = t18 * t54 - t32 * t55;
t4 = t14 * t55 + t24 * t54;
t3 = t14 * t54 - t24 * t55;
t2 = t12 * t55 + t23 * t54;
t1 = t12 * t54 - t23 * t55;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t65, -t62, 0, 0; t62, t65, 0, 0; 0, 0, 1, t93; 0, 0, 0, 1; t36, t35, t99, t92; t34, t33, -t98, t86; t45, t41, t67, t91; 0, 0, 0, 1; t14, -t13, t24, t85; t12, -t11, t23, t78; t18, -t17, t32, t84; 0, 0, 0, 1; t14 * t73 + t103, -t14 * t69 + t24 * t73, t13, t14 * pkin(3) + t13 * pkin(10) + t85; t12 * t73 + t104, -t12 * t69 + t23 * t73, t11, t12 * pkin(3) + t11 * pkin(10) + t78; t18 * t73 + t102, -t18 * t69 + t32 * t73, t17, t18 * pkin(3) + t17 * pkin(10) + t84; 0, 0, 0, 1; t4, -t3, t13, t81; t2, -t1, t11, t77; t6, -t5, t17, t79; 0, 0, 0, 1; t13 * t68 + t4 * t72, t13 * t72 - t4 * t68, t3, t4 * pkin(5) + t3 * pkin(12) + t81; t11 * t68 + t2 * t72, t11 * t72 - t2 * t68, t1, t2 * pkin(5) + t1 * pkin(12) + t77; t17 * t68 + t6 * t72, t17 * t72 - t6 * t68, t5, t6 * pkin(5) + t5 * pkin(12) + t79; 0, 0, 0, 1;];
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
