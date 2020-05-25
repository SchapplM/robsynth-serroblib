% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 17:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPRR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:58:59
% EndTime: 2018-11-23 17:59:00
% DurationCPUTime: 0.19s
% Computational Cost: add. (671->76), mult. (789->82), div. (0->0), fcn. (900->16), ass. (0->59)
t79 = pkin(9) - pkin(10);
t45 = sin(qJ(2));
t46 = sin(qJ(1));
t50 = cos(qJ(1));
t71 = pkin(6) - qJ(2);
t61 = cos(t71) / 0.2e1;
t70 = pkin(6) + qJ(2);
t65 = cos(t70);
t55 = t61 + t65 / 0.2e1;
t24 = t46 * t45 - t50 * t55;
t78 = t24 * pkin(9);
t26 = t50 * t45 + t46 * t55;
t77 = t26 * pkin(9);
t60 = sin(t70) / 0.2e1;
t64 = sin(t71);
t30 = t60 + t64 / 0.2e1;
t76 = t30 * pkin(9);
t75 = cos(qJ(3));
t41 = sin(pkin(6));
t74 = t46 * t41;
t73 = t50 * t41;
t72 = cos(pkin(6));
t69 = pkin(7) + 0;
t68 = t41 * t75;
t67 = t72 * pkin(8) + t69;
t66 = t50 * pkin(1) + pkin(8) * t74 + 0;
t32 = t61 - t65 / 0.2e1;
t63 = t32 * pkin(2) + t67;
t31 = t60 - t64 / 0.2e1;
t49 = cos(qJ(2));
t27 = -t46 * t31 + t50 * t49;
t62 = t27 * pkin(2) + t66;
t59 = t46 * pkin(1) - pkin(8) * t73 + 0;
t25 = t50 * t31 + t46 * t49;
t58 = t25 * pkin(2) + t59;
t44 = sin(qJ(3));
t22 = t32 * t44 - t72 * t75;
t23 = t32 * t75 + t72 * t44;
t57 = t23 * pkin(3) + t22 * qJ(4) + t63;
t15 = t27 * t44 - t46 * t68;
t16 = t27 * t75 + t44 * t74;
t56 = t16 * pkin(3) + t15 * qJ(4) + t62;
t13 = t25 * t44 + t50 * t68;
t14 = t25 * t75 - t44 * t73;
t54 = t14 * pkin(3) + t13 * qJ(4) + t58;
t53 = t23 * pkin(4) - t79 * t30 + t57;
t52 = t16 * pkin(4) + t79 * t26 + t56;
t51 = t14 * pkin(4) + t79 * t24 + t54;
t48 = cos(qJ(5));
t47 = cos(qJ(6));
t43 = sin(qJ(5));
t42 = sin(qJ(6));
t6 = t22 * t43 + t23 * t48;
t5 = -t22 * t48 + t23 * t43;
t4 = t15 * t43 + t16 * t48;
t3 = -t15 * t48 + t16 * t43;
t2 = t13 * t43 + t14 * t48;
t1 = -t13 * t48 + t14 * t43;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t50, -t46, 0, 0; t46, t50, 0, 0; 0, 0, 1, t69; 0, 0, 0, 1; t27, -t26, t74, t66; t25, -t24, -t73, t59; t32, t30, t72, t67; 0, 0, 0, 1; t16, -t15, t26, t62 + t77; t14, -t13, t24, t58 + t78; t23, -t22, -t30, t63 - t76; 0, 0, 0, 1; t16, t26, t15, t56 + t77; t14, t24, t13, t54 + t78; t23, -t30, t22, t57 - t76; 0, 0, 0, 1; t4, -t3, -t26, t52; t2, -t1, -t24, t51; t6, -t5, t30, t53; 0, 0, 0, 1; -t26 * t42 + t4 * t47, -t26 * t47 - t4 * t42, t3, t4 * pkin(5) + t3 * pkin(11) + t52; t2 * t47 - t24 * t42, -t2 * t42 - t24 * t47, t1, t2 * pkin(5) + t1 * pkin(11) + t51; t30 * t42 + t6 * t47, t30 * t47 - t6 * t42, t5, t6 * pkin(5) + t5 * pkin(11) + t53; 0, 0, 0, 1;];
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
