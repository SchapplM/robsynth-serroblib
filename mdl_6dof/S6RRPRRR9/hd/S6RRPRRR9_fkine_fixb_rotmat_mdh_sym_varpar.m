% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:26:31
% EndTime: 2018-11-23 17:26:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (541->89), mult. (472->103), div. (0->0), fcn. (534->18), ass. (0->55)
t42 = sin(pkin(12));
t72 = t42 * pkin(3);
t44 = cos(pkin(12));
t32 = t44 * pkin(3) + pkin(2);
t45 = cos(pkin(6));
t71 = t45 * t42;
t43 = sin(pkin(6));
t49 = sin(qJ(1));
t70 = t49 * t43;
t52 = cos(qJ(1));
t69 = t52 * t43;
t46 = -pkin(9) - qJ(3);
t68 = pkin(6) - qJ(2);
t67 = pkin(6) + qJ(2);
t66 = pkin(7) + 0;
t41 = pkin(12) + qJ(4);
t65 = t49 * pkin(1) + 0;
t64 = t42 * t70;
t63 = t45 * pkin(8) + t66;
t62 = t52 * pkin(1) + pkin(8) * t70 + 0;
t61 = cos(t67);
t60 = sin(t68);
t59 = cos(t68) / 0.2e1;
t58 = sin(t67) / 0.2e1;
t57 = -pkin(8) * t69 + t65;
t34 = cos(t41);
t20 = pkin(4) * t34 + t32;
t22 = t58 + t60 / 0.2e1;
t24 = t59 - t61 / 0.2e1;
t33 = sin(t41);
t25 = pkin(4) * t33 + t72;
t40 = -pkin(10) + t46;
t56 = t24 * t20 + t22 * t40 + t45 * t25 + t63;
t48 = sin(qJ(2));
t53 = t59 + t61 / 0.2e1;
t14 = t52 * t48 + t49 * t53;
t23 = t58 - t60 / 0.2e1;
t51 = cos(qJ(2));
t15 = -t49 * t23 + t52 * t51;
t55 = -t14 * t40 + t15 * t20 + t25 * t70 + t62;
t12 = t49 * t48 - t52 * t53;
t13 = t52 * t23 + t49 * t51;
t54 = t13 * t20 - t12 * t40 + (-pkin(8) - t25) * t69 + t65;
t50 = cos(qJ(6));
t47 = sin(qJ(6));
t35 = qJ(5) + t41;
t31 = cos(t35);
t30 = sin(t35);
t8 = t24 * t31 + t45 * t30;
t7 = t24 * t30 - t45 * t31;
t4 = t15 * t31 + t30 * t70;
t3 = t15 * t30 - t31 * t70;
t2 = t13 * t31 - t30 * t69;
t1 = t13 * t30 + t31 * t69;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t52, -t49, 0, 0; t49, t52, 0, 0; 0, 0, 1, t66; 0, 0, 0, 1; t15, -t14, t70, t62; t13, -t12, -t69, t57; t24, t22, t45, t63; 0, 0, 0, 1; t15 * t44 + t64, -t15 * t42 + t44 * t70, t14, t15 * pkin(2) + t14 * qJ(3) + t62; t13 * t44 - t42 * t69, -t13 * t42 - t44 * t69, t12, t13 * pkin(2) + t12 * qJ(3) + t57; t24 * t44 + t71, -t24 * t42 + t45 * t44, -t22, t24 * pkin(2) - t22 * qJ(3) + t63; 0, 0, 0, 1; t15 * t34 + t33 * t70, -t15 * t33 + t34 * t70, t14, pkin(3) * t64 - t14 * t46 + t15 * t32 + t62; t13 * t34 - t33 * t69, -t13 * t33 - t34 * t69, t12, -t12 * t46 + t13 * t32 + (-pkin(8) - t72) * t69 + t65; t24 * t34 + t45 * t33, -t24 * t33 + t45 * t34, -t22, pkin(3) * t71 + t22 * t46 + t24 * t32 + t63; 0, 0, 0, 1; t4, -t3, t14, t55; t2, -t1, t12, t54; t8, -t7, -t22, t56; 0, 0, 0, 1; t14 * t47 + t4 * t50, t14 * t50 - t4 * t47, t3, t4 * pkin(5) + t3 * pkin(11) + t55; t12 * t47 + t2 * t50, t12 * t50 - t2 * t47, t1, t2 * pkin(5) + t1 * pkin(11) + t54; -t22 * t47 + t8 * t50, -t22 * t50 - t8 * t47, t7, t8 * pkin(5) + t7 * pkin(11) + t56; 0, 0, 0, 1;];
T_ges = t5;
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
