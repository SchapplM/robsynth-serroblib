% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PPRRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:51:44
% EndTime: 2018-11-23 14:51:44
% DurationCPUTime: 0.31s
% Computational Cost: add. (1430->93), mult. (1432->117), div. (0->0), fcn. (1486->24), ass. (0->76)
t63 = sin(pkin(12));
t65 = sin(pkin(6));
t100 = t63 * t65;
t60 = pkin(6) - pkin(13);
t50 = cos(t60) / 0.2e1;
t59 = pkin(6) + pkin(13);
t53 = cos(t59);
t42 = t50 + t53 / 0.2e1;
t62 = sin(pkin(13));
t67 = cos(pkin(12));
t34 = -t63 * t42 - t67 * t62;
t64 = sin(pkin(7));
t68 = cos(pkin(7));
t24 = t68 * t100 - t34 * t64;
t49 = sin(t59) / 0.2e1;
t52 = sin(t60);
t40 = t49 + t52 / 0.2e1;
t69 = cos(pkin(6));
t36 = -t40 * t64 + t69 * t68;
t32 = t67 * t42 - t63 * t62;
t99 = t67 * t65;
t23 = -t32 * t64 - t68 * t99;
t71 = sin(qJ(4));
t105 = t23 * t71;
t104 = t24 * t71;
t102 = t36 * t71;
t97 = qJ(2) * t65;
t96 = pkin(7) - qJ(3);
t95 = pkin(7) + qJ(3);
t93 = qJ(1) + 0;
t92 = t67 * pkin(1) + t63 * t97 + 0;
t91 = t69 * qJ(2) + t93;
t90 = cos(t95);
t89 = sin(t96);
t88 = cos(t96) / 0.2e1;
t87 = sin(t95) / 0.2e1;
t86 = t63 * pkin(1) - t67 * t97 + 0;
t41 = t49 - t52 / 0.2e1;
t66 = cos(pkin(13));
t35 = -t63 * t41 + t67 * t66;
t85 = t35 * pkin(2) + t24 * pkin(8) + t92;
t43 = t50 - t53 / 0.2e1;
t84 = t43 * pkin(2) + t36 * pkin(8) + t91;
t83 = t88 + t90 / 0.2e1;
t82 = t87 + t89 / 0.2e1;
t72 = sin(qJ(3));
t80 = t65 * t82;
t13 = -t34 * t83 + t35 * t72 - t63 * t80;
t44 = t87 - t89 / 0.2e1;
t45 = t88 - t90 / 0.2e1;
t75 = cos(qJ(3));
t14 = t45 * t100 + t34 * t44 + t35 * t75;
t74 = cos(qJ(4));
t51 = t74 * pkin(4) + pkin(3);
t76 = -pkin(10) - pkin(9);
t81 = pkin(4) * t104 - t13 * t76 + t14 * t51 + t85;
t17 = -t40 * t83 + t43 * t72 - t69 * t82;
t18 = t40 * t44 + t43 * t75 + t69 * t45;
t79 = pkin(4) * t102 - t17 * t76 + t18 * t51 + t84;
t33 = t67 * t41 + t63 * t66;
t78 = t33 * pkin(2) + t23 * pkin(8) + t86;
t11 = -t32 * t83 + t33 * t72 + t67 * t80;
t12 = t32 * t44 + t33 * t75 - t45 * t99;
t77 = pkin(4) * t105 - t11 * t76 + t12 * t51 + t78;
t73 = cos(qJ(6));
t70 = sin(qJ(6));
t61 = qJ(4) + qJ(5);
t56 = cos(t61);
t55 = sin(t61);
t6 = t18 * t56 + t36 * t55;
t5 = t18 * t55 - t36 * t56;
t4 = t14 * t56 + t24 * t55;
t3 = t14 * t55 - t24 * t56;
t2 = t12 * t56 + t23 * t55;
t1 = t12 * t55 - t23 * t56;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t67, -t63, 0, 0; t63, t67, 0, 0; 0, 0, 1, t93; 0, 0, 0, 1; t35, t34, t100, t92; t33, t32, -t99, t86; t43, t40, t69, t91; 0, 0, 0, 1; t14, -t13, t24, t85; t12, -t11, t23, t78; t18, -t17, t36, t84; 0, 0, 0, 1; t14 * t74 + t104, -t14 * t71 + t24 * t74, t13, t14 * pkin(3) + t13 * pkin(9) + t85; t12 * t74 + t105, -t12 * t71 + t23 * t74, t11, t12 * pkin(3) + t11 * pkin(9) + t78; t18 * t74 + t102, -t18 * t71 + t36 * t74, t17, t18 * pkin(3) + t17 * pkin(9) + t84; 0, 0, 0, 1; t4, -t3, t13, t81; t2, -t1, t11, t77; t6, -t5, t17, t79; 0, 0, 0, 1; t13 * t70 + t4 * t73, t13 * t73 - t4 * t70, t3, t4 * pkin(5) + t3 * pkin(11) + t81; t11 * t70 + t2 * t73, t11 * t73 - t2 * t70, t1, t2 * pkin(5) + t1 * pkin(11) + t77; t17 * t70 + t6 * t73, t17 * t73 - t6 * t70, t5, t6 * pkin(5) + t5 * pkin(11) + t79; 0, 0, 0, 1;];
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
