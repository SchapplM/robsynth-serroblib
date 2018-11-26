% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:43:37
% EndTime: 2018-11-23 18:43:37
% DurationCPUTime: 0.29s
% Computational Cost: add. (1430->93), mult. (1432->116), div. (0->0), fcn. (1486->24), ass. (0->75)
t60 = pkin(6) - qJ(2);
t50 = cos(t60) / 0.2e1;
t59 = pkin(6) + qJ(2);
t53 = cos(t59);
t44 = t50 + t53 / 0.2e1;
t69 = sin(qJ(2));
t70 = sin(qJ(1));
t75 = cos(qJ(1));
t35 = -t70 * t44 - t75 * t69;
t62 = sin(pkin(7));
t64 = cos(pkin(7));
t63 = sin(pkin(6));
t98 = t70 * t63;
t24 = -t35 * t62 + t64 * t98;
t49 = sin(t59) / 0.2e1;
t52 = sin(t60);
t41 = t49 + t52 / 0.2e1;
t65 = cos(pkin(6));
t32 = -t41 * t62 + t65 * t64;
t33 = t75 * t44 - t70 * t69;
t97 = t75 * t63;
t23 = -t33 * t62 - t64 * t97;
t67 = sin(qJ(4));
t104 = t23 * t67;
t103 = t24 * t67;
t102 = t32 * t67;
t96 = pkin(7) - qJ(3);
t95 = pkin(7) + qJ(3);
t94 = pkin(8) + 0;
t92 = t65 * pkin(9) + t94;
t91 = t75 * pkin(1) + pkin(9) * t98 + 0;
t90 = cos(t95);
t89 = sin(t96);
t88 = cos(t96) / 0.2e1;
t87 = sin(t95) / 0.2e1;
t86 = t70 * pkin(1) - pkin(9) * t97 + 0;
t45 = t50 - t53 / 0.2e1;
t85 = t45 * pkin(2) + t32 * pkin(10) + t92;
t42 = t49 - t52 / 0.2e1;
t74 = cos(qJ(2));
t36 = -t70 * t42 + t75 * t74;
t84 = t36 * pkin(2) + t24 * pkin(10) + t91;
t83 = t88 + t90 / 0.2e1;
t82 = t87 + t89 / 0.2e1;
t68 = sin(qJ(3));
t17 = -t41 * t83 + t45 * t68 - t65 * t82;
t40 = t87 - t89 / 0.2e1;
t43 = t88 - t90 / 0.2e1;
t73 = cos(qJ(3));
t18 = t41 * t40 + t65 * t43 + t45 * t73;
t72 = cos(qJ(4));
t51 = t72 * pkin(4) + pkin(3);
t76 = -pkin(12) - pkin(11);
t81 = pkin(4) * t102 - t17 * t76 + t18 * t51 + t85;
t79 = t63 * t82;
t13 = -t35 * t83 + t36 * t68 - t70 * t79;
t14 = t35 * t40 + t36 * t73 + t43 * t98;
t80 = pkin(4) * t103 - t13 * t76 + t14 * t51 + t84;
t34 = t75 * t42 + t70 * t74;
t78 = t34 * pkin(2) + t23 * pkin(10) + t86;
t11 = -t33 * t83 + t34 * t68 + t75 * t79;
t12 = t33 * t40 + t34 * t73 - t43 * t97;
t77 = pkin(4) * t104 - t11 * t76 + t12 * t51 + t78;
t71 = cos(qJ(6));
t66 = sin(qJ(6));
t61 = qJ(4) + qJ(5);
t55 = cos(t61);
t54 = sin(t61);
t6 = t18 * t55 + t32 * t54;
t5 = t18 * t54 - t32 * t55;
t4 = t14 * t55 + t24 * t54;
t3 = t14 * t54 - t24 * t55;
t2 = t12 * t55 + t23 * t54;
t1 = t12 * t54 - t23 * t55;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t75, -t70, 0, 0; t70, t75, 0, 0; 0, 0, 1, t94; 0, 0, 0, 1; t36, t35, t98, t91; t34, t33, -t97, t86; t45, t41, t65, t92; 0, 0, 0, 1; t14, -t13, t24, t84; t12, -t11, t23, t78; t18, -t17, t32, t85; 0, 0, 0, 1; t14 * t72 + t103, -t14 * t67 + t24 * t72, t13, t14 * pkin(3) + t13 * pkin(11) + t84; t12 * t72 + t104, -t12 * t67 + t23 * t72, t11, t12 * pkin(3) + t11 * pkin(11) + t78; t18 * t72 + t102, -t18 * t67 + t32 * t72, t17, t18 * pkin(3) + t17 * pkin(11) + t85; 0, 0, 0, 1; t4, -t3, t13, t80; t2, -t1, t11, t77; t6, -t5, t17, t81; 0, 0, 0, 1; t13 * t66 + t4 * t71, t13 * t71 - t4 * t66, t3, t4 * pkin(5) + t3 * pkin(13) + t80; t11 * t66 + t2 * t71, t11 * t71 - t2 * t66, t1, t2 * pkin(5) + t1 * pkin(13) + t77; t17 * t66 + t6 * t71, t17 * t71 - t6 * t66, t5, t6 * pkin(5) + t5 * pkin(13) + t81; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
