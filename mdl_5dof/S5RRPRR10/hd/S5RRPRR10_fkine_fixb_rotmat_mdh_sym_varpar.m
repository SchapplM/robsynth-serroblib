% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:38
% EndTime: 2019-12-31 20:23:38
% DurationCPUTime: 0.18s
% Computational Cost: add. (213->58), mult. (485->79), div. (0->0), fcn. (668->12), ass. (0->46)
t32 = sin(pkin(10));
t34 = cos(pkin(10));
t38 = sin(qJ(2));
t42 = cos(qJ(2));
t21 = -t38 * t32 + t42 * t34;
t33 = sin(pkin(5));
t61 = t33 * t38;
t39 = sin(qJ(1));
t27 = t39 * t33;
t59 = t39 * t38;
t58 = t39 * t42;
t43 = cos(qJ(1));
t56 = t43 * t33;
t55 = t43 * t38;
t54 = t43 * t42;
t53 = pkin(6) + 0;
t35 = cos(pkin(5));
t52 = t35 * pkin(7) + t53;
t19 = t35 * t38 * pkin(2) + (-pkin(7) - qJ(3)) * t33;
t29 = t42 * pkin(2) + pkin(1);
t51 = t43 * t19 + t39 * t29 + 0;
t50 = t42 * t32 + t38 * t34;
t49 = -t39 * t19 + t43 * t29 + 0;
t48 = pkin(2) * t61 + t35 * qJ(3) + t52;
t47 = t21 * t35;
t7 = -t39 * t50 + t43 * t47;
t18 = t50 * t35;
t8 = t43 * t18 + t39 * t21;
t46 = t8 * pkin(3) - t7 * pkin(8) + t51;
t10 = -t39 * t18 + t43 * t21;
t9 = -t39 * t47 - t43 * t50;
t45 = t10 * pkin(3) - t9 * pkin(8) + t49;
t16 = t21 * t33;
t17 = t50 * t33;
t44 = t17 * pkin(3) - t16 * pkin(8) + t48;
t41 = cos(qJ(4));
t40 = cos(qJ(5));
t37 = sin(qJ(4));
t36 = sin(qJ(5));
t12 = t17 * t41 + t35 * t37;
t11 = t17 * t37 - t35 * t41;
t4 = t10 * t41 + t37 * t27;
t3 = t10 * t37 - t41 * t27;
t2 = -t37 * t56 + t8 * t41;
t1 = t8 * t37 + t41 * t56;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t43, -t39, 0, 0; t39, t43, 0, 0; 0, 0, 1, t53; 0, 0, 0, 1; -t35 * t59 + t54, -t35 * t58 - t55, t27, t43 * pkin(1) + pkin(7) * t27 + 0; t35 * t55 + t58, t35 * t54 - t59, -t56, t39 * pkin(1) - pkin(7) * t56 + 0; t61, t33 * t42, t35, t52; 0, 0, 0, 1; t10, t9, t27, t49; t8, t7, -t56, t51; t17, t16, t35, t48; 0, 0, 0, 1; t4, -t3, -t9, t45; t2, -t1, -t7, t46; t12, -t11, -t16, t44; 0, 0, 0, 1; -t9 * t36 + t4 * t40, -t4 * t36 - t9 * t40, t3, t4 * pkin(4) + t3 * pkin(9) + t45; t2 * t40 - t7 * t36, -t2 * t36 - t7 * t40, t1, t2 * pkin(4) + t1 * pkin(9) + t46; t12 * t40 - t16 * t36, -t12 * t36 - t16 * t40, t11, t12 * pkin(4) + t11 * pkin(9) + t44; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
