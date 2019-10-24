% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR4
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
% Datum: 2019-10-24 10:26
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:26:40
% EndTime: 2019-10-24 10:26:40
% DurationCPUTime: 0.18s
% Computational Cost: add. (213->58), mult. (485->80), div. (0->0), fcn. (668->12), ass. (0->44)
t32 = sin(pkin(10));
t35 = cos(pkin(10));
t40 = sin(qJ(2));
t43 = cos(qJ(2));
t21 = -t40 * t32 + t43 * t35;
t33 = sin(pkin(9));
t34 = sin(pkin(5));
t27 = t33 * t34;
t59 = t34 * t40;
t36 = cos(pkin(9));
t58 = t36 * t34;
t37 = cos(pkin(5));
t57 = t37 * t40;
t56 = t37 * t43;
t53 = qJ(1) + 0;
t19 = pkin(2) * t57 + (-pkin(6) - qJ(3)) * t34;
t29 = t43 * pkin(2) + pkin(1);
t52 = t36 * t19 + t33 * t29 + 0;
t51 = t37 * pkin(6) + t53;
t50 = t43 * t32 + t40 * t35;
t49 = -t33 * t19 + t36 * t29 + 0;
t48 = pkin(2) * t59 + t37 * qJ(3) + t51;
t47 = t21 * t37;
t7 = -t33 * t50 + t36 * t47;
t18 = t50 * t37;
t8 = t36 * t18 + t33 * t21;
t46 = t8 * pkin(3) - t7 * pkin(7) + t52;
t10 = -t33 * t18 + t36 * t21;
t9 = -t33 * t47 - t36 * t50;
t45 = t10 * pkin(3) - t9 * pkin(7) + t49;
t16 = t21 * t34;
t17 = t50 * t34;
t44 = t17 * pkin(3) - t16 * pkin(7) + t48;
t42 = cos(qJ(4));
t41 = cos(qJ(5));
t39 = sin(qJ(4));
t38 = sin(qJ(5));
t12 = t17 * t42 + t37 * t39;
t11 = t17 * t39 - t37 * t42;
t4 = t10 * t42 + t39 * t27;
t3 = t10 * t39 - t42 * t27;
t2 = -t39 * t58 + t8 * t42;
t1 = t8 * t39 + t42 * t58;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t36, -t33, 0, 0; t33, t36, 0, 0; 0, 0, 1, t53; 0, 0, 0, 1; -t33 * t57 + t36 * t43, -t33 * t56 - t36 * t40, t27, t36 * pkin(1) + pkin(6) * t27 + 0; t33 * t43 + t36 * t57, -t33 * t40 + t36 * t56, -t58, t33 * pkin(1) - pkin(6) * t58 + 0; t59, t34 * t43, t37, t51; 0, 0, 0, 1; t10, t9, t27, t49; t8, t7, -t58, t52; t17, t16, t37, t48; 0, 0, 0, 1; t4, -t3, -t9, t45; t2, -t1, -t7, t46; t12, -t11, -t16, t44; 0, 0, 0, 1; -t9 * t38 + t4 * t41, -t4 * t38 - t9 * t41, t3, t4 * pkin(4) + t3 * pkin(8) + t45; t2 * t41 - t7 * t38, -t2 * t38 - t7 * t41, t1, t2 * pkin(4) + t1 * pkin(8) + t46; t12 * t41 - t16 * t38, -t12 * t38 - t16 * t41, t11, t12 * pkin(4) + t11 * pkin(8) + t44; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
