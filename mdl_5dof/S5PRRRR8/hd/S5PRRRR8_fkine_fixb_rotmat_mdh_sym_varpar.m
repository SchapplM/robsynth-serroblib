% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:38
% EndTime: 2019-12-05 17:14:38
% DurationCPUTime: 0.15s
% Computational Cost: add. (188->62), mult. (307->80), div. (0->0), fcn. (422->12), ass. (0->43)
t30 = sin(pkin(10));
t31 = sin(pkin(5));
t57 = t30 * t31;
t36 = sin(qJ(2));
t56 = t31 * t36;
t38 = cos(qJ(3));
t55 = t31 * t38;
t39 = cos(qJ(2));
t54 = t31 * t39;
t32 = cos(pkin(10));
t53 = t32 * t31;
t33 = cos(pkin(5));
t35 = sin(qJ(3));
t52 = t33 * t35;
t51 = t33 * t36;
t50 = t33 * t39;
t49 = t30 * pkin(1) + 0;
t48 = t35 * t57;
t47 = qJ(1) + 0;
t46 = t32 * pkin(1) + pkin(6) * t57 + 0;
t45 = t33 * pkin(6) + t47;
t44 = -pkin(6) * t53 + t49;
t13 = t30 * t50 + t32 * t36;
t14 = -t30 * t51 + t32 * t39;
t23 = t38 * pkin(3) + pkin(2);
t40 = -pkin(8) - pkin(7);
t43 = pkin(3) * t48 - t13 * t40 + t14 * t23 + t46;
t42 = pkin(3) * t52 + t23 * t56 + t40 * t54 + t45;
t11 = t30 * t36 - t32 * t50;
t12 = t30 * t39 + t32 * t51;
t41 = t12 * t23 - t11 * t40 + (-pkin(3) * t35 - pkin(6)) * t53 + t49;
t37 = cos(qJ(5));
t34 = sin(qJ(5));
t29 = qJ(3) + qJ(4);
t25 = cos(t29);
t24 = sin(t29);
t8 = t33 * t24 + t25 * t56;
t7 = t24 * t56 - t33 * t25;
t4 = t14 * t25 + t24 * t57;
t3 = t14 * t24 - t25 * t57;
t2 = t12 * t25 - t24 * t53;
t1 = t12 * t24 + t25 * t53;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t32, -t30, 0, 0; t30, t32, 0, 0; 0, 0, 1, t47; 0, 0, 0, 1; t14, -t13, t57, t46; t12, -t11, -t53, t44; t56, t54, t33, t45; 0, 0, 0, 1; t14 * t38 + t48, -t14 * t35 + t30 * t55, t13, t14 * pkin(2) + t13 * pkin(7) + t46; t12 * t38 - t35 * t53, -t12 * t35 - t38 * t53, t11, t12 * pkin(2) + t11 * pkin(7) + t44; t36 * t55 + t52, t33 * t38 - t35 * t56, -t54, (pkin(2) * t36 - pkin(7) * t39) * t31 + t45; 0, 0, 0, 1; t4, -t3, t13, t43; t2, -t1, t11, t41; t8, -t7, -t54, t42; 0, 0, 0, 1; t13 * t34 + t4 * t37, t13 * t37 - t4 * t34, t3, t4 * pkin(4) + t3 * pkin(9) + t43; t11 * t34 + t2 * t37, t11 * t37 - t2 * t34, t1, t2 * pkin(4) + t1 * pkin(9) + t41; -t34 * t54 + t8 * t37, -t8 * t34 - t37 * t54, t7, t8 * pkin(4) + t7 * pkin(9) + t42; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
