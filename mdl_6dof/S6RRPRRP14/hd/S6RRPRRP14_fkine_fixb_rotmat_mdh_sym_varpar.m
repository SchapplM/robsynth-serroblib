% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRRP14_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:20:16
% EndTime: 2018-11-23 17:20:16
% DurationCPUTime: 0.16s
% Computational Cost: add. (550->70), mult. (613->69), div. (0->0), fcn. (686->14), ass. (0->53)
t41 = sin(pkin(6));
t46 = sin(qJ(1));
t36 = t46 * t41;
t50 = cos(qJ(1));
t73 = t50 * t41;
t72 = pkin(6) - qJ(2);
t71 = pkin(6) + qJ(2);
t70 = pkin(7) + 0;
t69 = pkin(8) * t73;
t68 = t46 * pkin(1) + 0;
t42 = cos(pkin(6));
t67 = t42 * pkin(8) + t70;
t66 = t50 * pkin(1) + pkin(8) * t36 + 0;
t65 = cos(t71);
t64 = sin(t72);
t63 = cos(t72) / 0.2e1;
t62 = sin(t71) / 0.2e1;
t61 = t63 + t65 / 0.2e1;
t45 = sin(qJ(2));
t21 = t46 * t45 - t50 * t61;
t49 = cos(qJ(2));
t57 = t62 - t64 / 0.2e1;
t22 = t46 * t49 + t50 * t57;
t60 = t22 * pkin(2) + t21 * qJ(3) + t68;
t30 = t62 + t64 / 0.2e1;
t31 = t63 - t65 / 0.2e1;
t59 = t31 * pkin(2) - t30 * qJ(3) + t67;
t23 = t50 * t45 + t46 * t61;
t24 = -t46 * t57 + t50 * t49;
t58 = t24 * pkin(2) + t23 * qJ(3) + t66;
t56 = t42 * pkin(3) + t31 * pkin(9) + t59;
t55 = pkin(3) * t36 + t24 * pkin(9) + t58;
t54 = t22 * pkin(9) + (-pkin(3) - pkin(8)) * t73 + t60;
t44 = sin(qJ(4));
t48 = cos(qJ(4));
t10 = t23 * t44 + t48 * t36;
t9 = -t23 * t48 + t44 * t36;
t53 = t10 * pkin(4) + t9 * pkin(10) + t55;
t19 = t30 * t48 + t42 * t44;
t20 = -t30 * t44 + t42 * t48;
t52 = t20 * pkin(4) + t19 * pkin(10) + t56;
t11 = t21 * t48 + t44 * t73;
t12 = t21 * t44 - t48 * t73;
t51 = t12 * pkin(4) - t11 * pkin(10) + t54;
t47 = cos(qJ(5));
t43 = sin(qJ(5));
t6 = t20 * t47 + t31 * t43;
t5 = t20 * t43 - t31 * t47;
t4 = t12 * t47 + t22 * t43;
t3 = t12 * t43 - t22 * t47;
t2 = t10 * t47 + t24 * t43;
t1 = t10 * t43 - t24 * t47;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t50, -t46, 0, 0; t46, t50, 0, 0; 0, 0, 1, t70; 0, 0, 0, 1; t24, -t23, t36, t66; t22, -t21, -t73, t68 - t69; t31, t30, t42, t67; 0, 0, 0, 1; t36, -t24, t23, t58; -t73, -t22, t21, t60 - t69; t42, -t31, -t30, t59; 0, 0, 0, 1; t10, -t9, t24, t55; t12, t11, t22, t54; t20, -t19, t31, t56; 0, 0, 0, 1; t2, -t1, t9, t53; t4, -t3, -t11, t51; t6, -t5, t19, t52; 0, 0, 0, 1; t2, t9, t1, t2 * pkin(5) + t1 * qJ(6) + t53; t4, -t11, t3, t4 * pkin(5) + t3 * qJ(6) + t51; t6, t19, t5, t6 * pkin(5) + t5 * qJ(6) + t52; 0, 0, 0, 1;];
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
