% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:46
% EndTime: 2019-12-31 20:16:46
% DurationCPUTime: 0.10s
% Computational Cost: add. (122->42), mult. (69->44), div. (0->0), fcn. (110->10), ass. (0->29)
t20 = sin(qJ(1));
t15 = qJ(2) + pkin(9);
t11 = qJ(4) + t15;
t6 = sin(t11);
t34 = t20 * t6;
t23 = cos(qJ(1));
t33 = t23 * t6;
t22 = cos(qJ(2));
t8 = t22 * pkin(2) + pkin(1);
t18 = sin(qJ(5));
t32 = t20 * t18;
t21 = cos(qJ(5));
t31 = t20 * t21;
t30 = t23 * t18;
t29 = t23 * t21;
t17 = -qJ(3) - pkin(6);
t16 = pkin(5) + 0;
t14 = -pkin(7) + t17;
t10 = cos(t15);
t3 = pkin(3) * t10 + t8;
t28 = t23 * t14 + t20 * t3 + 0;
t19 = sin(qJ(2));
t27 = t19 * pkin(2) + t16;
t9 = sin(t15);
t26 = pkin(3) * t9 + t27;
t7 = cos(t11);
t25 = pkin(4) * t7 + pkin(8) * t6;
t24 = -t20 * t14 + t23 * t3 + 0;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t20, 0, 0; t20, t23, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t23 * t22, -t23 * t19, t20, t23 * pkin(1) + t20 * pkin(6) + 0; t20 * t22, -t20 * t19, -t23, t20 * pkin(1) - t23 * pkin(6) + 0; t19, t22, 0, t16; 0, 0, 0, 1; t23 * t10, -t23 * t9, t20, -t20 * t17 + t23 * t8 + 0; t20 * t10, -t20 * t9, -t23, t23 * t17 + t20 * t8 + 0; t9, t10, 0, t27; 0, 0, 0, 1; t23 * t7, -t33, t20, t24; t20 * t7, -t34, -t23, t28; t6, t7, 0, t26; 0, 0, 0, 1; t7 * t29 + t32, -t7 * t30 + t31, t33, t25 * t23 + t24; t7 * t31 - t30, -t7 * t32 - t29, t34, t25 * t20 + t28; t6 * t21, -t6 * t18, -t7, t6 * pkin(4) - t7 * pkin(8) + t26; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
