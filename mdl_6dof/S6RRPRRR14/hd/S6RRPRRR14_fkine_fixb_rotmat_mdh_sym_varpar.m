% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRRR14_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:07:39
% EndTime: 2019-01-03 10:07:39
% DurationCPUTime: 0.30s
% Computational Cost: add. (901->82), mult. (2488->121), div. (0->0), fcn. (3324->18), ass. (0->74)
t59 = sin(pkin(6));
t71 = cos(qJ(2));
t100 = t59 * t71;
t58 = sin(pkin(7));
t62 = cos(pkin(7));
t63 = cos(pkin(6));
t44 = -t58 * t100 + t63 * t62;
t67 = sin(qJ(2));
t72 = cos(qJ(1));
t93 = t72 * t67;
t68 = sin(qJ(1));
t95 = t68 * t71;
t47 = -t63 * t95 - t93;
t97 = t68 * t59;
t39 = -t47 * t58 + t62 * t97;
t102 = t58 * t63;
t56 = sin(pkin(14));
t60 = cos(pkin(14));
t99 = t62 * t71;
t36 = t60 * t102 + (-t56 * t67 + t60 * t99) * t59;
t57 = sin(pkin(8));
t61 = cos(pkin(8));
t23 = -t36 * t57 + t44 * t61;
t92 = t72 * t71;
t96 = t68 * t67;
t48 = -t63 * t96 + t92;
t83 = t47 * t62 + t58 * t97;
t28 = -t48 * t56 + t83 * t60;
t19 = -t28 * t57 + t39 * t61;
t46 = t63 * t93 + t95;
t45 = t63 * t92 - t96;
t94 = t72 * t59;
t84 = t45 * t62 - t58 * t94;
t26 = -t46 * t56 + t84 * t60;
t38 = -t45 * t58 - t62 * t94;
t18 = -t26 * t57 + t38 * t61;
t110 = cos(qJ(4));
t101 = t59 * t67;
t91 = pkin(9) + 0;
t88 = t57 * t110;
t87 = t61 * t110;
t86 = t63 * pkin(10) + t91;
t85 = t72 * pkin(1) + pkin(10) * t97 + 0;
t82 = t68 * pkin(1) - pkin(10) * t94 + 0;
t81 = t48 * pkin(2) + t39 * qJ(3) + t85;
t80 = pkin(2) * t101 + t44 * qJ(3) + t86;
t29 = t48 * t60 + t83 * t56;
t79 = t29 * pkin(3) + t19 * pkin(11) + t81;
t37 = t60 * t101 + (t59 * t99 + t102) * t56;
t78 = t37 * pkin(3) + t23 * pkin(11) + t80;
t77 = t46 * pkin(2) + t38 * qJ(3) + t82;
t66 = sin(qJ(4));
t11 = -t28 * t87 + t29 * t66 - t39 * t88;
t12 = t29 * t110 + (t28 * t61 + t39 * t57) * t66;
t76 = t12 * pkin(4) + t11 * pkin(12) + t79;
t14 = -t36 * t87 + t37 * t66 - t44 * t88;
t15 = t37 * t110 + (t36 * t61 + t44 * t57) * t66;
t75 = t15 * pkin(4) + t14 * pkin(12) + t78;
t27 = t46 * t60 + t84 * t56;
t74 = t27 * pkin(3) + t18 * pkin(11) + t77;
t10 = t27 * t110 + (t26 * t61 + t38 * t57) * t66;
t9 = -t26 * t87 + t27 * t66 - t38 * t88;
t73 = t10 * pkin(4) + t9 * pkin(12) + t74;
t70 = cos(qJ(5));
t69 = cos(qJ(6));
t65 = sin(qJ(5));
t64 = sin(qJ(6));
t6 = t15 * t70 + t23 * t65;
t5 = t15 * t65 - t23 * t70;
t4 = t12 * t70 + t19 * t65;
t3 = t12 * t65 - t19 * t70;
t2 = t10 * t70 + t18 * t65;
t1 = t10 * t65 - t18 * t70;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t72, -t68, 0, 0; t68, t72, 0, 0; 0, 0, 1, t91; 0, 0, 0, 1; t48, t47, t97, t85; t46, t45, -t94, t82; t101, t100, t63, t86; 0, 0, 0, 1; t29, t28, t39, t81; t27, t26, t38, t77; t37, t36, t44, t80; 0, 0, 0, 1; t12, -t11, t19, t79; t10, -t9, t18, t74; t15, -t14, t23, t78; 0, 0, 0, 1; t4, -t3, t11, t76; t2, -t1, t9, t73; t6, -t5, t14, t75; 0, 0, 0, 1; t11 * t64 + t4 * t69, t11 * t69 - t4 * t64, t3, t4 * pkin(5) + t3 * pkin(13) + t76; t2 * t69 + t9 * t64, -t2 * t64 + t9 * t69, t1, t2 * pkin(5) + t1 * pkin(13) + t73; t14 * t64 + t6 * t69, t14 * t69 - t6 * t64, t5, t6 * pkin(5) + t5 * pkin(13) + t75; 0, 0, 0, 1;];
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
