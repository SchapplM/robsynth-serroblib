% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:42
% EndTime: 2019-12-31 21:52:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (111->43), mult. (92->44), div. (0->0), fcn. (138->8), ass. (0->33)
t17 = qJ(2) + qJ(3);
t13 = sin(t17);
t19 = sin(qJ(4));
t36 = t13 * t19;
t21 = sin(qJ(1));
t35 = t21 * t19;
t22 = cos(qJ(4));
t34 = t21 * t22;
t24 = cos(qJ(1));
t33 = t24 * t19;
t32 = t24 * t22;
t16 = pkin(5) + 0;
t23 = cos(qJ(2));
t11 = t23 * pkin(2) + pkin(1);
t31 = t24 * t11 + 0;
t20 = sin(qJ(2));
t30 = t20 * pkin(2) + t16;
t25 = -pkin(7) - pkin(6);
t29 = t21 * t11 + t24 * t25 + 0;
t14 = cos(t17);
t28 = pkin(3) * t14 + pkin(8) * t13;
t10 = t22 * pkin(4) + pkin(3);
t18 = -qJ(5) - pkin(8);
t27 = t10 * t14 - t13 * t18;
t26 = -t21 * t25 + t31;
t9 = t24 * t13;
t8 = t21 * t13;
t7 = t13 * t22;
t4 = t14 * t32 + t35;
t3 = -t14 * t33 + t34;
t2 = t14 * t34 - t33;
t1 = -t14 * t35 - t32;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t21, 0, 0; t21, t24, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t24 * t23, -t24 * t20, t21, t24 * pkin(1) + t21 * pkin(6) + 0; t21 * t23, -t21 * t20, -t24, t21 * pkin(1) - t24 * pkin(6) + 0; t20, t23, 0, t16; 0, 0, 0, 1; t24 * t14, -t9, t21, t26; t21 * t14, -t8, -t24, t29; t13, t14, 0, t30; 0, 0, 0, 1; t4, t3, t9, t24 * t28 + t26; t2, t1, t8, t21 * t28 + t29; t7, -t36, -t14, t13 * pkin(3) - t14 * pkin(8) + t30; 0, 0, 0, 1; t4, t3, t9, t27 * t24 + (pkin(4) * t19 - t25) * t21 + t31; t2, t1, t8, -pkin(4) * t33 + t21 * t27 + t29; t7, -t36, -t14, t13 * t10 + t14 * t18 + t30; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
