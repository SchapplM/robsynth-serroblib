% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-10-24 10:22
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRPPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:22:49
% EndTime: 2019-10-24 10:22:49
% DurationCPUTime: 0.13s
% Computational Cost: add. (107->44), mult. (80->40), div. (0->0), fcn. (121->8), ass. (0->29)
t18 = cos(pkin(7));
t16 = qJ(2) + pkin(8);
t12 = sin(t16);
t30 = qJ(4) * t12;
t13 = cos(t16);
t8 = t18 * t13;
t37 = pkin(3) * t8 + t18 * t30;
t17 = sin(pkin(7));
t36 = t17 * t12;
t7 = t17 * t13;
t20 = sin(qJ(5));
t35 = t17 * t20;
t22 = cos(qJ(5));
t34 = t17 * t22;
t33 = t18 * t12;
t32 = t18 * t20;
t31 = t18 * t22;
t23 = cos(qJ(2));
t11 = t23 * pkin(2) + pkin(1);
t29 = t18 * t11 + 0;
t15 = qJ(1) + 0;
t19 = -qJ(3) - pkin(5);
t28 = t17 * t11 + t18 * t19 + 0;
t21 = sin(qJ(2));
t27 = t21 * pkin(2) + t15;
t26 = pkin(3) * t7 + t17 * t30 + t28;
t25 = -t17 * t19 + t29;
t24 = t12 * pkin(3) - t13 * qJ(4) + t27;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t17, 0, 0; t17, t18, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t18 * t23, -t18 * t21, t17, t18 * pkin(1) + t17 * pkin(5) + 0; t17 * t23, -t17 * t21, -t18, t17 * pkin(1) - t18 * pkin(5) + 0; t21, t23, 0, t15; 0, 0, 0, 1; t8, -t33, t17, t25; t7, -t36, -t18, t28; t12, t13, 0, t27; 0, 0, 0, 1; t17, -t8, t33, t25 + t37; -t18, -t7, t36, t26; 0, -t12, -t13, t24; 0, 0, 0, 1; t12 * t32 + t34, t12 * t31 - t35, t8, pkin(6) * t8 + (pkin(4) - t19) * t17 + t29 + t37; t12 * t35 - t31, t12 * t34 + t32, t7, -t18 * pkin(4) + pkin(6) * t7 + t26; -t13 * t20, -t13 * t22, t12, t12 * pkin(6) + t24; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
