% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:39
% EndTime: 2019-12-31 19:56:39
% DurationCPUTime: 0.12s
% Computational Cost: add. (111->43), mult. (92->44), div. (0->0), fcn. (138->8), ass. (0->33)
t16 = qJ(2) + pkin(8);
t13 = sin(t16);
t20 = sin(qJ(4));
t36 = t13 * t20;
t22 = sin(qJ(1));
t35 = t22 * t20;
t23 = cos(qJ(4));
t34 = t22 * t23;
t25 = cos(qJ(1));
t33 = t25 * t20;
t32 = t25 * t23;
t17 = pkin(5) + 0;
t24 = cos(qJ(2));
t12 = t24 * pkin(2) + pkin(1);
t31 = t25 * t12 + 0;
t21 = sin(qJ(2));
t30 = t21 * pkin(2) + t17;
t19 = -qJ(3) - pkin(6);
t29 = t22 * t12 + t25 * t19 + 0;
t14 = cos(t16);
t28 = pkin(3) * t14 + pkin(7) * t13;
t11 = t23 * pkin(4) + pkin(3);
t18 = -qJ(5) - pkin(7);
t27 = t11 * t14 - t13 * t18;
t26 = -t22 * t19 + t31;
t9 = t25 * t13;
t8 = t22 * t13;
t7 = t13 * t23;
t4 = t14 * t32 + t35;
t3 = -t14 * t33 + t34;
t2 = t14 * t34 - t33;
t1 = -t14 * t35 - t32;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, -t22, 0, 0; t22, t25, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; t25 * t24, -t25 * t21, t22, t25 * pkin(1) + t22 * pkin(6) + 0; t22 * t24, -t22 * t21, -t25, t22 * pkin(1) - t25 * pkin(6) + 0; t21, t24, 0, t17; 0, 0, 0, 1; t25 * t14, -t9, t22, t26; t22 * t14, -t8, -t25, t29; t13, t14, 0, t30; 0, 0, 0, 1; t4, t3, t9, t28 * t25 + t26; t2, t1, t8, t28 * t22 + t29; t7, -t36, -t14, t13 * pkin(3) - t14 * pkin(7) + t30; 0, 0, 0, 1; t4, t3, t9, t27 * t25 + (pkin(4) * t20 - t19) * t22 + t31; t2, t1, t8, -pkin(4) * t33 + t27 * t22 + t29; t7, -t36, -t14, t13 * t11 + t14 * t18 + t30; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
