% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRRP10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:00
% EndTime: 2019-12-31 22:08:00
% DurationCPUTime: 0.13s
% Computational Cost: add. (167->54), mult. (371->64), div. (0->0), fcn. (511->10), ass. (0->48)
t28 = sin(pkin(5));
t33 = sin(qJ(2));
t57 = t28 * t33;
t37 = cos(qJ(2));
t56 = t28 * t37;
t34 = sin(qJ(1));
t55 = t34 * t28;
t54 = t34 * t33;
t53 = t34 * t37;
t38 = cos(qJ(1));
t52 = t38 * t28;
t51 = t38 * t33;
t50 = t38 * t37;
t49 = pkin(6) + 0;
t31 = sin(qJ(4));
t48 = pkin(4) * t31 + pkin(8);
t29 = cos(pkin(5));
t47 = t29 * pkin(7) + t49;
t46 = t38 * pkin(1) + pkin(7) * t55 + 0;
t45 = pkin(2) * t57 + t47;
t18 = -t29 * t54 + t50;
t44 = t18 * pkin(2) + t46;
t43 = t34 * pkin(1) - pkin(7) * t52 + 0;
t16 = t29 * t51 + t53;
t42 = t16 * pkin(2) + t43;
t17 = t29 * t53 + t51;
t41 = t17 * pkin(8) + t44;
t40 = -pkin(8) * t56 + t45;
t15 = -t29 * t50 + t54;
t39 = t15 * pkin(8) + t42;
t36 = cos(qJ(3));
t35 = cos(qJ(4));
t32 = sin(qJ(3));
t30 = -qJ(5) - pkin(9);
t24 = t35 * pkin(4) + pkin(3);
t14 = t29 * t32 + t36 * t57;
t13 = -t29 * t36 + t32 * t57;
t10 = t18 * t36 + t32 * t55;
t9 = t18 * t32 - t36 * t55;
t8 = t16 * t36 - t32 * t52;
t7 = t16 * t32 + t36 * t52;
t6 = t14 * t35 - t31 * t56;
t5 = -t14 * t31 - t35 * t56;
t4 = t10 * t35 + t17 * t31;
t3 = -t10 * t31 + t17 * t35;
t2 = t15 * t31 + t8 * t35;
t1 = t15 * t35 - t8 * t31;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t38, -t34, 0, 0; t34, t38, 0, 0; 0, 0, 1, t49; 0, 0, 0, 1; t18, -t17, t55, t46; t16, -t15, -t52, t43; t57, t56, t29, t47; 0, 0, 0, 1; t10, -t9, t17, t41; t8, -t7, t15, t39; t14, -t13, -t56, t40; 0, 0, 0, 1; t4, t3, t9, t10 * pkin(3) + t9 * pkin(9) + t41; t2, t1, t7, t8 * pkin(3) + t7 * pkin(9) + t39; t6, t5, t13, t14 * pkin(3) + t13 * pkin(9) + t40; 0, 0, 0, 1; t4, t3, t9, t10 * t24 + t48 * t17 - t9 * t30 + t44; t2, t1, t7, t48 * t15 + t8 * t24 - t7 * t30 + t42; t6, t5, t13, -t13 * t30 + t14 * t24 - t48 * t56 + t45; 0, 0, 0, 1;];
T_ges = t11;
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
