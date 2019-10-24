% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-10-24 10:28
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRPRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:27:56
% EndTime: 2019-10-24 10:27:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (88->51), mult. (110->52), div. (0->0), fcn. (156->8), ass. (0->29)
t16 = sin(pkin(8));
t19 = sin(qJ(2));
t31 = qJ(3) * t19;
t21 = cos(qJ(2));
t5 = t16 * t21;
t36 = pkin(2) * t5 + t16 * t31;
t35 = t16 * t19;
t17 = cos(pkin(8));
t34 = t17 * t19;
t6 = t17 * t21;
t18 = sin(qJ(4));
t33 = t18 * t19;
t20 = cos(qJ(4));
t32 = t19 * t20;
t30 = t16 * pkin(1) + 0;
t14 = qJ(1) + 0;
t29 = t17 * pkin(1) + t16 * pkin(5) + 0;
t28 = t19 * pkin(2) + t14;
t27 = t30 + t36;
t26 = -t17 * pkin(5) + t30;
t25 = pkin(2) * t6 + t17 * t31 + t29;
t22 = -pkin(7) - pkin(6);
t24 = pkin(4) * t33 - t21 * t22;
t23 = -t21 * qJ(3) + t28;
t15 = qJ(4) + qJ(5);
t9 = cos(t15);
t8 = sin(t15);
t7 = t20 * pkin(4) + pkin(3);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t16, 0, 0; t16, t17, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; t6, -t34, t16, t29; t5, -t35, -t17, t26; t19, t21, 0, t14; 0, 0, 0, 1; t16, -t6, t34, t25; -t17, -t5, t35, t26 + t36; 0, -t19, -t21, t23; 0, 0, 0, 1; t16 * t20 + t17 * t33, -t16 * t18 + t17 * t32, t6, t16 * pkin(3) + pkin(6) * t6 + t25; t16 * t33 - t17 * t20, t16 * t32 + t17 * t18, t5, pkin(6) * t5 + (-pkin(3) - pkin(5)) * t17 + t27; -t21 * t18, -t21 * t20, t19, t19 * pkin(6) + t23; 0, 0, 0, 1; t16 * t9 + t8 * t34, -t16 * t8 + t9 * t34, t6, t16 * t7 + t24 * t17 + t25; -t17 * t9 + t8 * t35, t17 * t8 + t9 * t35, t5, (-pkin(5) - t7) * t17 + t24 * t16 + t27; -t21 * t8, -t21 * t9, t19, -t19 * t22 + (-pkin(4) * t18 - qJ(3)) * t21 + t28; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
