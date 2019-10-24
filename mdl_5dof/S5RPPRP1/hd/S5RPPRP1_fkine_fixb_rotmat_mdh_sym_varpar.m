% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-10-24 10:39
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:39:39
% EndTime: 2019-10-24 10:39:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (123->43), mult. (83->37), div. (0->0), fcn. (125->8), ass. (0->32)
t15 = qJ(1) + pkin(7);
t11 = sin(t15);
t16 = sin(pkin(8));
t35 = t11 * t16;
t19 = sin(qJ(4));
t34 = t11 * t19;
t12 = cos(t15);
t33 = t12 * t19;
t32 = t16 * t19;
t17 = cos(pkin(8));
t31 = t17 * t19;
t21 = cos(qJ(4));
t30 = t17 * t21;
t29 = pkin(5) + 0;
t22 = cos(qJ(1));
t28 = t22 * pkin(1) + 0;
t13 = qJ(2) + t29;
t20 = sin(qJ(1));
t27 = -t20 * pkin(1) + 0;
t26 = t12 * pkin(2) + t11 * qJ(3) + t28;
t25 = pkin(3) * t17 + pkin(6) * t16;
t24 = t12 * qJ(3) + t27;
t10 = t21 * pkin(4) + pkin(3);
t18 = -qJ(5) - pkin(6);
t23 = t10 * t17 - t16 * t18;
t9 = t16 * t21;
t5 = t12 * t16;
t4 = t12 * t30 + t34;
t3 = t11 * t21 - t12 * t31;
t2 = -t11 * t30 + t33;
t1 = t11 * t31 + t12 * t21;
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t29; -t20, -t22, 0, 0; t22, -t20, 0, 0; 0, 0, 0, 1; 0, 0, 1, t13; -t11, -t12, 0, t27; t12, -t11, 0, t28; 0, 0, 0, 1; t16, t17, 0, t13; -t11 * t17, t35, t12, -t11 * pkin(2) + t24; t12 * t17, -t5, t11, t26; 0, 0, 0, 1; t9, -t32, -t17, t16 * pkin(3) - t17 * pkin(6) + t13; t2, t1, -t35, (-pkin(2) - t25) * t11 + t24; t4, t3, t5, t25 * t12 + t26; 0, 0, 0, 1; t9, -t32, -t17, t16 * t10 + t17 * t18 + t13; t2, t1, -t35, pkin(4) * t33 + (-pkin(2) - t23) * t11 + t24; t4, t3, t5, pkin(4) * t34 + t23 * t12 + t26; 0, 0, 0, 1;];
T_ges = t6;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
