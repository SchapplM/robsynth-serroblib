% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR9
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:19:54
% EndTime: 2019-12-31 20:19:54
% DurationCPUTime: 0.11s
% Computational Cost: add. (121->48), mult. (92->54), div. (0->0), fcn. (138->10), ass. (0->34)
t15 = qJ(4) + qJ(5);
t10 = sin(t15);
t19 = sin(qJ(1));
t37 = t19 * t10;
t11 = cos(t15);
t36 = t19 * t11;
t17 = sin(qJ(4));
t35 = t19 * t17;
t20 = cos(qJ(4));
t34 = t19 * t20;
t22 = cos(qJ(1));
t33 = t22 * t10;
t32 = t22 * t11;
t31 = t22 * t17;
t30 = t22 * t20;
t14 = pkin(5) + 0;
t21 = cos(qJ(2));
t7 = t21 * pkin(2) + pkin(1);
t29 = t22 * t7 + 0;
t16 = -qJ(3) - pkin(6);
t28 = t22 * t16 + t19 * t7 + 0;
t18 = sin(qJ(2));
t27 = t18 * pkin(2) + t14;
t13 = qJ(2) + pkin(9);
t8 = sin(t13);
t9 = cos(t13);
t26 = pkin(3) * t9 + pkin(7) * t8;
t23 = -pkin(8) - pkin(7);
t6 = t20 * pkin(4) + pkin(3);
t25 = -t23 * t8 + t6 * t9;
t24 = -t19 * t16 + t29;
t4 = t22 * t8;
t3 = t19 * t8;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t19, 0, 0; t19, t22, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; t22 * t21, -t22 * t18, t19, t22 * pkin(1) + t19 * pkin(6) + 0; t19 * t21, -t19 * t18, -t22, t19 * pkin(1) - t22 * pkin(6) + 0; t18, t21, 0, t14; 0, 0, 0, 1; t22 * t9, -t4, t19, t24; t19 * t9, -t3, -t22, t28; t8, t9, 0, t27; 0, 0, 0, 1; t9 * t30 + t35, -t9 * t31 + t34, t4, t26 * t22 + t24; t9 * t34 - t31, -t9 * t35 - t30, t3, t26 * t19 + t28; t8 * t20, -t8 * t17, -t9, t8 * pkin(3) - t9 * pkin(7) + t27; 0, 0, 0, 1; t9 * t32 + t37, -t9 * t33 + t36, t4, t25 * t22 + (pkin(4) * t17 - t16) * t19 + t29; t9 * t36 - t33, -t9 * t37 - t32, t3, -pkin(4) * t31 + t25 * t19 + t28; t8 * t11, -t8 * t10, -t9, t9 * t23 + t8 * t6 + t27; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
