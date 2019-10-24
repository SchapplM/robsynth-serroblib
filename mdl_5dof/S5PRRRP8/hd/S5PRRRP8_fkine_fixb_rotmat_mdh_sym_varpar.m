% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP8
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
% Datum: 2019-10-24 10:35
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRP8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:35:15
% EndTime: 2019-10-24 10:35:16
% DurationCPUTime: 0.15s
% Computational Cost: add. (182->51), mult. (417->62), div. (0->0), fcn. (575->10), ass. (0->44)
t32 = sin(pkin(9));
t33 = sin(pkin(5));
t58 = t32 * t33;
t38 = sin(qJ(2));
t57 = t33 * t38;
t40 = cos(qJ(3));
t56 = t33 * t40;
t41 = cos(qJ(2));
t55 = t33 * t41;
t34 = cos(pkin(9));
t54 = t34 * t33;
t35 = cos(pkin(5));
t53 = t35 * t38;
t52 = t35 * t41;
t51 = qJ(1) + 0;
t50 = t34 * pkin(1) + pkin(6) * t58 + 0;
t49 = t35 * pkin(6) + t51;
t48 = t32 * pkin(1) - pkin(6) * t54 + 0;
t20 = t32 * t52 + t34 * t38;
t21 = -t32 * t53 + t34 * t41;
t47 = t21 * pkin(2) + t20 * pkin(7) + t50;
t46 = pkin(2) * t57 - pkin(7) * t55 + t49;
t18 = t32 * t38 - t34 * t52;
t19 = t32 * t41 + t34 * t53;
t45 = t19 * pkin(2) + t18 * pkin(7) + t48;
t37 = sin(qJ(3));
t10 = t21 * t40 + t37 * t58;
t9 = t21 * t37 - t32 * t56;
t44 = t10 * pkin(3) + t9 * pkin(8) + t47;
t22 = -t35 * t40 + t37 * t57;
t23 = t35 * t37 + t38 * t56;
t43 = t23 * pkin(3) + t22 * pkin(8) + t46;
t7 = t19 * t37 + t40 * t54;
t8 = t19 * t40 - t37 * t54;
t42 = t8 * pkin(3) + t7 * pkin(8) + t45;
t39 = cos(qJ(4));
t36 = sin(qJ(4));
t12 = t23 * t39 - t36 * t55;
t11 = t23 * t36 + t39 * t55;
t4 = t10 * t39 + t20 * t36;
t3 = t10 * t36 - t20 * t39;
t2 = t18 * t36 + t8 * t39;
t1 = -t18 * t39 + t8 * t36;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t32, 0, 0; t32, t34, 0, 0; 0, 0, 1, t51; 0, 0, 0, 1; t21, -t20, t58, t50; t19, -t18, -t54, t48; t57, t55, t35, t49; 0, 0, 0, 1; t10, -t9, t20, t47; t8, -t7, t18, t45; t23, -t22, -t55, t46; 0, 0, 0, 1; t4, -t3, t9, t44; t2, -t1, t7, t42; t12, -t11, t22, t43; 0, 0, 0, 1; t4, t9, t3, t4 * pkin(4) + t3 * qJ(5) + t44; t2, t7, t1, t2 * pkin(4) + t1 * qJ(5) + t42; t12, t22, t11, t12 * pkin(4) + t11 * qJ(5) + t43; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
