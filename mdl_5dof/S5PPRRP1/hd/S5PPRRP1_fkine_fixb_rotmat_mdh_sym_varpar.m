% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRP1
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
% Datum: 2019-10-24 10:19
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PPRRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:19:15
% EndTime: 2019-10-24 10:19:15
% DurationCPUTime: 0.13s
% Computational Cost: add. (111->43), mult. (92->44), div. (0->0), fcn. (138->8), ass. (0->33)
t17 = pkin(8) + qJ(3);
t13 = sin(t17);
t24 = sin(qJ(4));
t36 = t13 * t24;
t19 = sin(pkin(7));
t35 = t19 * t24;
t25 = cos(qJ(4));
t34 = t19 * t25;
t21 = cos(pkin(7));
t33 = t21 * t24;
t32 = t21 * t25;
t20 = cos(pkin(8));
t11 = t20 * pkin(2) + pkin(1);
t31 = t21 * t11 + 0;
t16 = qJ(1) + 0;
t23 = -pkin(5) - qJ(2);
t30 = t19 * t11 + t21 * t23 + 0;
t18 = sin(pkin(8));
t29 = t18 * pkin(2) + t16;
t14 = cos(t17);
t28 = pkin(3) * t14 + pkin(6) * t13;
t12 = t25 * pkin(4) + pkin(3);
t22 = -qJ(5) - pkin(6);
t27 = t12 * t14 - t13 * t22;
t26 = -t19 * t23 + t31;
t9 = t13 * t25;
t8 = t21 * t13;
t7 = t19 * t13;
t4 = t14 * t32 + t35;
t3 = -t14 * t33 + t34;
t2 = t14 * t34 - t33;
t1 = -t14 * t35 - t32;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t19, 0, 0; t19, t21, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t21 * t20, -t21 * t18, t19, t21 * pkin(1) + t19 * qJ(2) + 0; t19 * t20, -t19 * t18, -t21, t19 * pkin(1) - t21 * qJ(2) + 0; t18, t20, 0, t16; 0, 0, 0, 1; t21 * t14, -t8, t19, t26; t19 * t14, -t7, -t21, t30; t13, t14, 0, t29; 0, 0, 0, 1; t4, t3, t8, t28 * t21 + t26; t2, t1, t7, t28 * t19 + t30; t9, -t36, -t14, t13 * pkin(3) - t14 * pkin(6) + t29; 0, 0, 0, 1; t4, t3, t8, t27 * t21 + (pkin(4) * t24 - t23) * t19 + t31; t2, t1, t7, -pkin(4) * t33 + t27 * t19 + t30; t9, -t36, -t14, t13 * t12 + t14 * t22 + t29; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
