% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PPPRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:17:37
% EndTime: 2019-10-24 10:17:38
% DurationCPUTime: 0.14s
% Computational Cost: add. (122->50), mult. (250->61), div. (0->0), fcn. (348->10), ass. (0->37)
t27 = sin(pkin(9));
t28 = sin(pkin(8));
t51 = t28 * t27;
t30 = cos(pkin(9));
t50 = t28 * t30;
t29 = sin(pkin(7));
t49 = t29 * t28;
t31 = cos(pkin(8));
t48 = t29 * t31;
t32 = cos(pkin(7));
t47 = t32 * t28;
t46 = t32 * t31;
t45 = qJ(3) * t28;
t26 = qJ(1) + 0;
t44 = t32 * pkin(1) + t29 * qJ(2) + 0;
t43 = t29 * pkin(1) - t32 * qJ(2) + 0;
t42 = pkin(2) * t46 + t32 * t45 + t44;
t41 = t28 * pkin(2) - t31 * qJ(3) + t26;
t40 = pkin(2) * t48 + t29 * t45 + t43;
t39 = pkin(3) * t50 + pkin(5) * t51 + t41;
t10 = t29 * t27 + t30 * t46;
t9 = t27 * t46 - t29 * t30;
t38 = t10 * pkin(3) + t9 * pkin(5) + t42;
t7 = t27 * t48 + t32 * t30;
t8 = -t32 * t27 + t30 * t48;
t37 = t8 * pkin(3) + t7 * pkin(5) + t40;
t36 = cos(qJ(4));
t35 = cos(qJ(5));
t34 = sin(qJ(4));
t33 = sin(qJ(5));
t12 = -t31 * t34 + t36 * t50;
t11 = t31 * t36 + t34 * t50;
t4 = t10 * t36 + t34 * t47;
t3 = t10 * t34 - t36 * t47;
t2 = t34 * t49 + t8 * t36;
t1 = t8 * t34 - t36 * t49;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t32, -t29, 0, 0; t29, t32, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t46, -t47, t29, t44; t48, -t49, -t32, t43; t28, t31, 0, t26; 0, 0, 0, 1; t10, -t9, t47, t42; t8, -t7, t49, t40; t50, -t51, -t31, t41; 0, 0, 0, 1; t4, -t3, t9, t38; t2, -t1, t7, t37; t12, -t11, t51, t39; 0, 0, 0, 1; t9 * t33 + t4 * t35, -t4 * t33 + t9 * t35, t3, t4 * pkin(4) + t3 * pkin(6) + t38; t2 * t35 + t7 * t33, -t2 * t33 + t7 * t35, t1, t2 * pkin(4) + t1 * pkin(6) + t37; t12 * t35 + t33 * t51, -t12 * t33 + t35 * t51, t11, t12 * pkin(4) + t11 * pkin(6) + t39; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
