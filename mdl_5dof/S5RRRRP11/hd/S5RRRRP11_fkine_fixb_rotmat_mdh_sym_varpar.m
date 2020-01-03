% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP11
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRRP11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:17
% EndTime: 2019-12-31 22:14:17
% DurationCPUTime: 0.15s
% Computational Cost: add. (182->51), mult. (417->59), div. (0->0), fcn. (575->10), ass. (0->45)
t32 = sin(pkin(5));
t36 = sin(qJ(2));
t59 = t32 * t36;
t40 = cos(qJ(2));
t58 = t32 * t40;
t37 = sin(qJ(1));
t57 = t37 * t32;
t56 = t37 * t36;
t55 = t37 * t40;
t41 = cos(qJ(1));
t54 = t41 * t32;
t53 = t41 * t36;
t52 = t41 * t40;
t51 = pkin(6) + 0;
t33 = cos(pkin(5));
t50 = t33 * pkin(7) + t51;
t49 = t41 * pkin(1) + pkin(7) * t57 + 0;
t48 = t37 * pkin(1) - pkin(7) * t54 + 0;
t22 = t33 * t55 + t53;
t23 = -t33 * t56 + t52;
t47 = t23 * pkin(2) + t22 * pkin(8) + t49;
t46 = pkin(2) * t59 - pkin(8) * t58 + t50;
t20 = -t33 * t52 + t56;
t21 = t33 * t53 + t55;
t45 = t21 * pkin(2) + t20 * pkin(8) + t48;
t35 = sin(qJ(3));
t39 = cos(qJ(3));
t11 = t23 * t35 - t39 * t57;
t12 = t23 * t39 + t35 * t57;
t44 = t12 * pkin(3) + t11 * pkin(9) + t47;
t18 = -t33 * t39 + t35 * t59;
t19 = t33 * t35 + t39 * t59;
t43 = t19 * pkin(3) + t18 * pkin(9) + t46;
t10 = t21 * t39 - t35 * t54;
t9 = t21 * t35 + t39 * t54;
t42 = t10 * pkin(3) + t9 * pkin(9) + t45;
t38 = cos(qJ(4));
t34 = sin(qJ(4));
t8 = t19 * t38 - t34 * t58;
t7 = t19 * t34 + t38 * t58;
t4 = t12 * t38 + t22 * t34;
t3 = t12 * t34 - t22 * t38;
t2 = t10 * t38 + t20 * t34;
t1 = t10 * t34 - t20 * t38;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t41, -t37, 0, 0; t37, t41, 0, 0; 0, 0, 1, t51; 0, 0, 0, 1; t23, -t22, t57, t49; t21, -t20, -t54, t48; t59, t58, t33, t50; 0, 0, 0, 1; t12, -t11, t22, t47; t10, -t9, t20, t45; t19, -t18, -t58, t46; 0, 0, 0, 1; t4, -t3, t11, t44; t2, -t1, t9, t42; t8, -t7, t18, t43; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(4) + t3 * qJ(5) + t44; t2, t9, t1, t2 * pkin(4) + t1 * qJ(5) + t42; t8, t18, t7, t8 * pkin(4) + t7 * qJ(5) + t43; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
