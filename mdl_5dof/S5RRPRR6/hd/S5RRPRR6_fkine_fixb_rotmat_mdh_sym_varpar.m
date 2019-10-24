% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-10-24 10:49
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:49:47
% EndTime: 2019-10-24 10:49:47
% DurationCPUTime: 0.12s
% Computational Cost: add. (133->48), mult. (83->47), div. (0->0), fcn. (125->10), ass. (0->31)
t14 = sin(pkin(9));
t13 = qJ(1) + qJ(2);
t7 = sin(t13);
t34 = t7 * t14;
t15 = cos(pkin(9));
t33 = t7 * t15;
t16 = sin(qJ(4));
t32 = t7 * t16;
t9 = cos(t13);
t31 = t9 * t15;
t30 = t9 * t16;
t29 = t15 * t16;
t18 = cos(qJ(4));
t28 = t15 * t18;
t27 = pkin(5) + 0;
t19 = cos(qJ(1));
t26 = t19 * pkin(1) + 0;
t10 = pkin(6) + t27;
t17 = sin(qJ(1));
t25 = -t17 * pkin(1) + 0;
t24 = t9 * pkin(2) + t7 * qJ(3) + t26;
t23 = pkin(3) * t15 + pkin(7) * t14;
t20 = -pkin(8) - pkin(7);
t5 = t18 * pkin(4) + pkin(3);
t22 = -t14 * t20 + t15 * t5;
t21 = t9 * qJ(3) + t25;
t12 = qJ(4) + qJ(5);
t8 = cos(t12);
t6 = sin(t12);
t1 = t9 * t14;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t27; -t17, -t19, 0, 0; t19, -t17, 0, 0; 0, 0, 0, 1; 0, 0, 1, t10; -t7, -t9, 0, t25; t9, -t7, 0, t26; 0, 0, 0, 1; t14, t15, 0, t10; -t33, t34, t9, -t7 * pkin(2) + t21; t31, -t1, t7, t24; 0, 0, 0, 1; t14 * t18, -t14 * t16, -t15, t14 * pkin(3) - t15 * pkin(7) + t10; -t7 * t28 + t30, t9 * t18 + t7 * t29, -t34, (-pkin(2) - t23) * t7 + t21; t9 * t28 + t32, t7 * t18 - t9 * t29, t1, t23 * t9 + t24; 0, 0, 0, 1; t14 * t8, -t14 * t6, -t15, t14 * t5 + t15 * t20 + t10; -t8 * t33 + t9 * t6, t6 * t33 + t9 * t8, -t34, pkin(4) * t30 + (-pkin(2) - t22) * t7 + t21; t8 * t31 + t7 * t6, -t6 * t31 + t7 * t8, t1, pkin(4) * t32 + t22 * t9 + t24; 0, 0, 0, 1;];
T_ges = t2;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
