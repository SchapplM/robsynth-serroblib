% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-10-24 10:46
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:46:34
% EndTime: 2019-10-24 10:46:34
% DurationCPUTime: 0.16s
% Computational Cost: add. (124->60), mult. (117->64), div. (0->0), fcn. (168->10), ass. (0->32)
t23 = -pkin(7) - pkin(6);
t21 = cos(qJ(3));
t6 = t21 * pkin(3) + pkin(2);
t17 = sin(pkin(9));
t20 = sin(qJ(1));
t35 = t20 * t17;
t18 = cos(pkin(9));
t34 = t20 * t18;
t19 = sin(qJ(3));
t33 = t20 * t19;
t32 = t20 * t21;
t22 = cos(qJ(1));
t31 = t22 * t18;
t30 = t22 * t19;
t29 = t22 * t21;
t16 = qJ(3) + qJ(4);
t14 = pkin(5) + 0;
t28 = t22 * qJ(2) + 0;
t27 = t22 * pkin(1) + t20 * qJ(2) + 0;
t26 = pkin(2) * t18 + pkin(6) * t17;
t8 = cos(t16);
t1 = pkin(4) * t8 + t6;
t15 = -pkin(8) + t23;
t25 = t1 * t18 - t15 * t17;
t24 = -t17 * t23 + t18 * t6;
t11 = qJ(5) + t16;
t7 = sin(t16);
t5 = cos(t11);
t4 = sin(t11);
t3 = t22 * t17;
t2 = t19 * pkin(3) + pkin(4) * t7;
t9 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t14; -t20, -t22, 0, 0; t22, -t20, 0, 0; 0, 0, 0, 1; t17, t18, 0, t14; -t34, t35, t22, -t20 * pkin(1) + t28; t31, -t3, t20, t27; 0, 0, 0, 1; t17 * t21, -t17 * t19, -t18, t17 * pkin(2) - t18 * pkin(6) + t14; -t18 * t32 + t30, t18 * t33 + t29, -t35, (-pkin(1) - t26) * t20 + t28; t18 * t29 + t33, -t18 * t30 + t32, t3, t26 * t22 + t27; 0, 0, 0, 1; t17 * t8, -t17 * t7, -t18, t17 * t6 + t18 * t23 + t14; t22 * t7 - t8 * t34, t22 * t8 + t7 * t34, -t35, pkin(3) * t30 + (-pkin(1) - t24) * t20 + t28; t20 * t7 + t8 * t31, t20 * t8 - t7 * t31, t3, pkin(3) * t33 + t24 * t22 + t27; 0, 0, 0, 1; t17 * t5, -t17 * t4, -t18, t17 * t1 + t18 * t15 + t14; t22 * t4 - t5 * t34, t22 * t5 + t4 * t34, -t35, t22 * t2 + (-pkin(1) - t25) * t20 + t28; t20 * t4 + t5 * t31, t20 * t5 - t4 * t31, t3, t20 * t2 + t25 * t22 + t27; 0, 0, 0, 1;];
T_ges = t9;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
