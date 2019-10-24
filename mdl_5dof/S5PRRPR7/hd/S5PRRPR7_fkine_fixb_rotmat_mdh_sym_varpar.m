% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-10-24 10:32
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRPR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:32:10
% EndTime: 2019-10-24 10:32:10
% DurationCPUTime: 0.14s
% Computational Cost: add. (203->57), mult. (472->74), div. (0->0), fcn. (651->12), ass. (0->46)
t33 = sin(pkin(9));
t34 = sin(pkin(5));
t60 = t33 * t34;
t40 = sin(qJ(2));
t59 = t34 * t40;
t42 = cos(qJ(3));
t58 = t34 * t42;
t43 = cos(qJ(2));
t57 = t34 * t43;
t36 = cos(pkin(9));
t56 = t36 * t34;
t37 = cos(pkin(5));
t55 = t37 * t40;
t54 = t37 * t43;
t53 = qJ(1) + 0;
t52 = t36 * pkin(1) + pkin(6) * t60 + 0;
t51 = t37 * pkin(6) + t53;
t50 = t33 * pkin(1) - pkin(6) * t56 + 0;
t20 = t33 * t54 + t36 * t40;
t21 = -t33 * t55 + t36 * t43;
t49 = t21 * pkin(2) + t20 * pkin(7) + t52;
t48 = pkin(2) * t59 - pkin(7) * t57 + t51;
t18 = t33 * t40 - t36 * t54;
t19 = t33 * t43 + t36 * t55;
t47 = t19 * pkin(2) + t18 * pkin(7) + t50;
t39 = sin(qJ(3));
t11 = t21 * t39 - t33 * t58;
t12 = t21 * t42 + t39 * t60;
t46 = t12 * pkin(3) + t11 * qJ(4) + t49;
t22 = -t37 * t42 + t39 * t59;
t23 = t37 * t39 + t40 * t58;
t45 = t23 * pkin(3) + t22 * qJ(4) + t48;
t10 = t19 * t42 - t39 * t56;
t9 = t19 * t39 + t42 * t56;
t44 = t10 * pkin(3) + t9 * qJ(4) + t47;
t41 = cos(qJ(5));
t38 = sin(qJ(5));
t35 = cos(pkin(10));
t32 = sin(pkin(10));
t8 = t23 * t35 - t32 * t57;
t7 = t23 * t32 + t35 * t57;
t4 = t12 * t35 + t20 * t32;
t3 = t12 * t32 - t20 * t35;
t2 = t10 * t35 + t18 * t32;
t1 = t10 * t32 - t18 * t35;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t36, -t33, 0, 0; t33, t36, 0, 0; 0, 0, 1, t53; 0, 0, 0, 1; t21, -t20, t60, t52; t19, -t18, -t56, t50; t59, t57, t37, t51; 0, 0, 0, 1; t12, -t11, t20, t49; t10, -t9, t18, t47; t23, -t22, -t57, t48; 0, 0, 0, 1; t4, -t3, t11, t46; t2, -t1, t9, t44; t8, -t7, t22, t45; 0, 0, 0, 1; t11 * t38 + t4 * t41, t11 * t41 - t4 * t38, t3, t4 * pkin(4) + t3 * pkin(8) + t46; t2 * t41 + t9 * t38, -t2 * t38 + t9 * t41, t1, t2 * pkin(4) + t1 * pkin(8) + t44; t22 * t38 + t8 * t41, t22 * t41 - t8 * t38, t7, t8 * pkin(4) + t7 * pkin(8) + t45; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
