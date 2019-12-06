% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:53:58
% EndTime: 2019-12-05 16:53:58
% DurationCPUTime: 0.14s
% Computational Cost: add. (167->54), mult. (371->67), div. (0->0), fcn. (511->10), ass. (0->47)
t28 = sin(pkin(9));
t29 = sin(pkin(5));
t56 = t28 * t29;
t35 = sin(qJ(2));
t55 = t29 * t35;
t37 = cos(qJ(3));
t54 = t29 * t37;
t38 = cos(qJ(2));
t53 = t29 * t38;
t30 = cos(pkin(9));
t52 = t30 * t29;
t31 = cos(pkin(5));
t51 = t31 * t35;
t50 = t31 * t38;
t49 = qJ(1) + 0;
t33 = sin(qJ(4));
t48 = pkin(4) * t33 + pkin(7);
t47 = t30 * pkin(1) + pkin(6) * t56 + 0;
t46 = t31 * pkin(6) + t49;
t16 = -t28 * t51 + t30 * t38;
t45 = t16 * pkin(2) + t47;
t44 = pkin(2) * t55 + t46;
t43 = t28 * pkin(1) - pkin(6) * t52 + 0;
t14 = t28 * t38 + t30 * t51;
t42 = t14 * pkin(2) + t43;
t15 = t28 * t50 + t30 * t35;
t41 = t15 * pkin(7) + t45;
t40 = -pkin(7) * t53 + t44;
t13 = t28 * t35 - t30 * t50;
t39 = t13 * pkin(7) + t42;
t36 = cos(qJ(4));
t34 = sin(qJ(3));
t32 = -qJ(5) - pkin(8);
t24 = t36 * pkin(4) + pkin(3);
t18 = t31 * t34 + t35 * t54;
t17 = -t31 * t37 + t34 * t55;
t10 = t18 * t36 - t33 * t53;
t9 = -t18 * t33 - t36 * t53;
t8 = t16 * t37 + t34 * t56;
t7 = t16 * t34 - t28 * t54;
t6 = t14 * t37 - t34 * t52;
t5 = t14 * t34 + t37 * t52;
t4 = t15 * t33 + t8 * t36;
t3 = t15 * t36 - t8 * t33;
t2 = t13 * t33 + t6 * t36;
t1 = t13 * t36 - t6 * t33;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t28, 0, 0; t28, t30, 0, 0; 0, 0, 1, t49; 0, 0, 0, 1; t16, -t15, t56, t47; t14, -t13, -t52, t43; t55, t53, t31, t46; 0, 0, 0, 1; t8, -t7, t15, t41; t6, -t5, t13, t39; t18, -t17, -t53, t40; 0, 0, 0, 1; t4, t3, t7, t8 * pkin(3) + t7 * pkin(8) + t41; t2, t1, t5, t6 * pkin(3) + t5 * pkin(8) + t39; t10, t9, t17, t18 * pkin(3) + t17 * pkin(8) + t40; 0, 0, 0, 1; t4, t3, t7, t48 * t15 + t8 * t24 - t7 * t32 + t45; t2, t1, t5, t48 * t13 + t6 * t24 - t5 * t32 + t42; t10, t9, t17, -t17 * t32 + t18 * t24 - t48 * t53 + t44; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
