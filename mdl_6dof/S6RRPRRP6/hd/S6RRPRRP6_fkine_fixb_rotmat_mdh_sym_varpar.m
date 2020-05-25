% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:14:08
% EndTime: 2018-11-23 17:14:08
% DurationCPUTime: 0.24s
% Computational Cost: add. (771->78), mult. (589->86), div. (0->0), fcn. (643->20), ass. (0->64)
t55 = pkin(6) - qJ(2);
t43 = cos(t55) / 0.2e1;
t54 = pkin(6) + qJ(2);
t51 = cos(t54);
t87 = t43 - t51 / 0.2e1;
t42 = sin(t54) / 0.2e1;
t49 = sin(t55);
t32 = t42 - t49 / 0.2e1;
t56 = sin(pkin(6));
t62 = sin(qJ(1));
t44 = t62 * t56;
t66 = cos(qJ(1));
t84 = t66 * t56;
t83 = pkin(7) + 0;
t53 = qJ(2) + pkin(11);
t58 = pkin(8) + qJ(3);
t25 = t32 * pkin(2) - t56 * t58;
t65 = cos(qJ(2));
t46 = t65 * pkin(2) + pkin(1);
t82 = t66 * t25 + t62 * t46 + 0;
t81 = pkin(6) - t53;
t80 = pkin(6) + t53;
t79 = -t62 * t25 + t66 * t46 + 0;
t78 = cos(t80);
t77 = sin(t81);
t57 = cos(pkin(6));
t76 = t87 * pkin(2) + t57 * t58 + t83;
t75 = cos(t81) / 0.2e1;
t74 = sin(t80) / 0.2e1;
t47 = sin(t53);
t67 = t78 / 0.2e1 + t75;
t17 = t62 * t47 - t66 * t67;
t29 = t74 - t77 / 0.2e1;
t50 = cos(t53);
t18 = t66 * t29 + t62 * t50;
t73 = t18 * pkin(3) + t17 * pkin(9) + t82;
t19 = t66 * t47 + t62 * t67;
t20 = -t62 * t29 + t66 * t50;
t72 = t20 * pkin(3) + t19 * pkin(9) + t79;
t30 = t77 / 0.2e1 + t74;
t31 = t75 - t78 / 0.2e1;
t71 = t31 * pkin(3) - t30 * pkin(9) + t76;
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t10 = t18 * t64 - t60 * t84;
t9 = t18 * t60 + t64 * t84;
t70 = t10 * pkin(4) + t9 * pkin(10) + t73;
t11 = t20 * t60 - t64 * t44;
t12 = t20 * t64 + t60 * t44;
t69 = t12 * pkin(4) + t11 * pkin(10) + t72;
t22 = t31 * t60 - t57 * t64;
t23 = t31 * t64 + t57 * t60;
t68 = t23 * pkin(4) + t22 * pkin(10) + t71;
t63 = cos(qJ(5));
t61 = sin(qJ(2));
t59 = sin(qJ(5));
t33 = t43 + t51 / 0.2e1;
t6 = t23 * t63 - t30 * t59;
t5 = t23 * t59 + t30 * t63;
t4 = t12 * t63 + t19 * t59;
t3 = t12 * t59 - t19 * t63;
t2 = t10 * t63 + t17 * t59;
t1 = t10 * t59 - t17 * t63;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t66, -t62, 0, 0; t62, t66, 0, 0; 0, 0, 1, t83; 0, 0, 0, 1; -t62 * t32 + t66 * t65, -t62 * t33 - t66 * t61, t44, t66 * pkin(1) + pkin(8) * t44 + 0; t66 * t32 + t62 * t65, t66 * t33 - t62 * t61, -t84, t62 * pkin(1) - pkin(8) * t84 + 0; t87, t42 + t49 / 0.2e1, t57, t57 * pkin(8) + t83; 0, 0, 0, 1; t20, -t19, t44, t79; t18, -t17, -t84, t82; t31, t30, t57, t76; 0, 0, 0, 1; t12, -t11, t19, t72; t10, -t9, t17, t73; t23, -t22, -t30, t71; 0, 0, 0, 1; t4, -t3, t11, t69; t2, -t1, t9, t70; t6, -t5, t22, t68; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(5) + t3 * qJ(6) + t69; t2, t9, t1, t2 * pkin(5) + t1 * qJ(6) + t70; t6, t22, t5, t6 * pkin(5) + t5 * qJ(6) + t68; 0, 0, 0, 1;];
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
