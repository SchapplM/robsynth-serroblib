% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRPRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:17
% EndTime: 2019-12-05 15:56:18
% DurationCPUTime: 0.14s
% Computational Cost: add. (188->62), mult. (307->79), div. (0->0), fcn. (422->12), ass. (0->42)
t31 = sin(pkin(9));
t32 = sin(pkin(5));
t56 = t31 * t32;
t38 = sin(qJ(2));
t55 = t32 * t38;
t40 = cos(qJ(2));
t54 = t32 * t40;
t34 = cos(pkin(9));
t53 = t34 * t32;
t30 = sin(pkin(10));
t35 = cos(pkin(5));
t52 = t35 * t30;
t51 = t35 * t38;
t50 = t35 * t40;
t49 = t31 * pkin(1) + 0;
t48 = t30 * t56;
t47 = qJ(1) + 0;
t46 = t34 * pkin(1) + pkin(6) * t56 + 0;
t45 = t35 * pkin(6) + t47;
t44 = -pkin(6) * t53 + t49;
t13 = t31 * t50 + t34 * t38;
t14 = -t31 * t51 + t34 * t40;
t33 = cos(pkin(10));
t23 = t33 * pkin(3) + pkin(2);
t36 = -pkin(7) - qJ(3);
t43 = pkin(3) * t48 - t13 * t36 + t14 * t23 + t46;
t42 = pkin(3) * t52 + t23 * t55 + t36 * t54 + t45;
t11 = t31 * t38 - t34 * t50;
t12 = t31 * t40 + t34 * t51;
t41 = t12 * t23 - t11 * t36 + (-pkin(3) * t30 - pkin(6)) * t53 + t49;
t39 = cos(qJ(5));
t37 = sin(qJ(5));
t29 = pkin(10) + qJ(4);
t25 = cos(t29);
t24 = sin(t29);
t8 = t35 * t24 + t25 * t55;
t7 = t24 * t55 - t35 * t25;
t4 = t14 * t25 + t24 * t56;
t3 = t14 * t24 - t25 * t56;
t2 = t12 * t25 - t24 * t53;
t1 = t12 * t24 + t25 * t53;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t31, 0, 0; t31, t34, 0, 0; 0, 0, 1, t47; 0, 0, 0, 1; t14, -t13, t56, t46; t12, -t11, -t53, t44; t55, t54, t35, t45; 0, 0, 0, 1; t14 * t33 + t48, -t14 * t30 + t33 * t56, t13, t14 * pkin(2) + t13 * qJ(3) + t46; t12 * t33 - t30 * t53, -t12 * t30 - t33 * t53, t11, t12 * pkin(2) + t11 * qJ(3) + t44; t33 * t55 + t52, -t30 * t55 + t35 * t33, -t54, (pkin(2) * t38 - qJ(3) * t40) * t32 + t45; 0, 0, 0, 1; t4, -t3, t13, t43; t2, -t1, t11, t41; t8, -t7, -t54, t42; 0, 0, 0, 1; t13 * t37 + t4 * t39, t13 * t39 - t4 * t37, t3, t4 * pkin(4) + t3 * pkin(8) + t43; t11 * t37 + t2 * t39, t11 * t39 - t2 * t37, t1, t2 * pkin(4) + t1 * pkin(8) + t41; -t37 * t54 + t8 * t39, -t8 * t37 - t39 * t54, t7, t8 * pkin(4) + t7 * pkin(8) + t42; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
