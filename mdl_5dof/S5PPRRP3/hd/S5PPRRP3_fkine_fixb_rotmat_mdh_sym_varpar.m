% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-10-24 10:20
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PPRRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:20:00
% EndTime: 2019-10-24 10:20:01
% DurationCPUTime: 0.11s
% Computational Cost: add. (109->44), mult. (216->48), div. (0->0), fcn. (301->8), ass. (0->38)
t28 = sin(pkin(8));
t33 = sin(qJ(3));
t21 = t28 * t33;
t35 = cos(qJ(3));
t52 = t28 * t35;
t29 = sin(pkin(7));
t51 = t29 * t28;
t30 = cos(pkin(8));
t50 = t29 * t30;
t49 = t29 * t33;
t48 = t29 * t35;
t31 = cos(pkin(7));
t47 = t31 * t28;
t46 = t31 * t30;
t45 = t31 * t33;
t44 = t31 * t35;
t27 = qJ(1) + 0;
t43 = t31 * pkin(1) + t29 * qJ(2) + 0;
t42 = t29 * pkin(1) - t31 * qJ(2) + 0;
t41 = pkin(2) * t46 + pkin(5) * t47 + t43;
t40 = t28 * pkin(2) - t30 * pkin(5) + t27;
t39 = pkin(2) * t50 + pkin(5) * t51 + t42;
t38 = pkin(3) * t52 + pkin(6) * t21 + t40;
t10 = t30 * t44 + t49;
t9 = t30 * t45 - t48;
t37 = t10 * pkin(3) + t9 * pkin(6) + t41;
t7 = t30 * t49 + t44;
t8 = t30 * t48 - t45;
t36 = t8 * pkin(3) + t7 * pkin(6) + t39;
t34 = cos(qJ(4));
t32 = sin(qJ(4));
t12 = -t30 * t32 + t34 * t52;
t11 = t30 * t34 + t32 * t52;
t4 = t10 * t34 + t32 * t47;
t3 = t10 * t32 - t34 * t47;
t2 = t32 * t51 + t8 * t34;
t1 = t8 * t32 - t34 * t51;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t31, -t29, 0, 0; t29, t31, 0, 0; 0, 0, 1, t27; 0, 0, 0, 1; t46, -t47, t29, t43; t50, -t51, -t31, t42; t28, t30, 0, t27; 0, 0, 0, 1; t10, -t9, t47, t41; t8, -t7, t51, t39; t52, -t21, -t30, t40; 0, 0, 0, 1; t4, -t3, t9, t37; t2, -t1, t7, t36; t12, -t11, t21, t38; 0, 0, 0, 1; t4, t9, t3, t4 * pkin(4) + t3 * qJ(5) + t37; t2, t7, t1, t2 * pkin(4) + t1 * qJ(5) + t36; t12, t21, t11, t12 * pkin(4) + t11 * qJ(5) + t38; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
