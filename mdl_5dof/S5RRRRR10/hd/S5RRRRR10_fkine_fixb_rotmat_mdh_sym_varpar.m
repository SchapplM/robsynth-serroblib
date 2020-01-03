% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRRR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:31:50
% EndTime: 2019-12-31 22:31:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (188->62), mult. (307->77), div. (0->0), fcn. (422->12), ass. (0->44)
t30 = sin(pkin(5));
t34 = sin(qJ(2));
t58 = t30 * t34;
t38 = cos(qJ(2));
t57 = t30 * t38;
t31 = cos(pkin(5));
t33 = sin(qJ(3));
t56 = t31 * t33;
t35 = sin(qJ(1));
t55 = t35 * t30;
t54 = t35 * t34;
t53 = t35 * t38;
t39 = cos(qJ(1));
t52 = t39 * t30;
t51 = t39 * t34;
t50 = t39 * t38;
t49 = pkin(6) + 0;
t48 = t35 * pkin(1) + 0;
t47 = t33 * t55;
t46 = t31 * pkin(7) + t49;
t45 = t39 * pkin(1) + pkin(7) * t55 + 0;
t44 = -pkin(7) * t52 + t48;
t37 = cos(qJ(3));
t23 = t37 * pkin(3) + pkin(2);
t40 = -pkin(9) - pkin(8);
t43 = pkin(3) * t56 + t23 * t58 + t40 * t57 + t46;
t13 = t31 * t53 + t51;
t14 = -t31 * t54 + t50;
t42 = pkin(3) * t47 - t13 * t40 + t14 * t23 + t45;
t11 = -t31 * t50 + t54;
t12 = t31 * t51 + t53;
t41 = t12 * t23 - t11 * t40 + (-pkin(3) * t33 - pkin(7)) * t52 + t48;
t36 = cos(qJ(5));
t32 = sin(qJ(5));
t29 = qJ(3) + qJ(4);
t25 = cos(t29);
t24 = sin(t29);
t8 = t31 * t24 + t25 * t58;
t7 = t24 * t58 - t31 * t25;
t4 = t14 * t25 + t24 * t55;
t3 = t14 * t24 - t25 * t55;
t2 = t12 * t25 - t24 * t52;
t1 = t12 * t24 + t25 * t52;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t39, -t35, 0, 0; t35, t39, 0, 0; 0, 0, 1, t49; 0, 0, 0, 1; t14, -t13, t55, t45; t12, -t11, -t52, t44; t58, t57, t31, t46; 0, 0, 0, 1; t14 * t37 + t47, -t14 * t33 + t37 * t55, t13, t14 * pkin(2) + t13 * pkin(8) + t45; t12 * t37 - t33 * t52, -t12 * t33 - t37 * t52, t11, t12 * pkin(2) + t11 * pkin(8) + t44; t37 * t58 + t56, t31 * t37 - t33 * t58, -t57, (pkin(2) * t34 - pkin(8) * t38) * t30 + t46; 0, 0, 0, 1; t4, -t3, t13, t42; t2, -t1, t11, t41; t8, -t7, -t57, t43; 0, 0, 0, 1; t13 * t32 + t4 * t36, t13 * t36 - t4 * t32, t3, t4 * pkin(4) + t3 * pkin(10) + t42; t11 * t32 + t2 * t36, t11 * t36 - t2 * t32, t1, t2 * pkin(4) + t1 * pkin(10) + t41; -t32 * t57 + t8 * t36, -t8 * t32 - t36 * t57, t7, t8 * pkin(4) + t7 * pkin(10) + t43; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
