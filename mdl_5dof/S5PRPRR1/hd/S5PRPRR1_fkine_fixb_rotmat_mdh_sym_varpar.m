% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-10-24 10:25
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:25:31
% EndTime: 2019-10-24 10:25:31
% DurationCPUTime: 0.09s
% Computational Cost: add. (110->36), mult. (41->30), div. (0->0), fcn. (73->10), ass. (0->23)
t20 = cos(pkin(9));
t4 = t20 * pkin(3) + pkin(2);
t22 = -pkin(6) - qJ(3);
t16 = pkin(9) + qJ(4);
t19 = sin(pkin(8));
t26 = t19 * pkin(1) + 0;
t21 = cos(pkin(8));
t25 = t21 * pkin(1) + 0;
t24 = qJ(1) + 0;
t10 = pkin(5) + t24;
t18 = sin(pkin(9));
t23 = t18 * pkin(3) + t10;
t17 = pkin(8) + qJ(2);
t15 = -pkin(7) + t22;
t9 = qJ(5) + t16;
t8 = cos(t17);
t7 = cos(t16);
t6 = sin(t17);
t5 = sin(t16);
t3 = cos(t9);
t2 = sin(t9);
t1 = pkin(4) * t7 + t4;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t19, 0, 0; t19, t21, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t8, -t6, 0, t25; t6, t8, 0, t26; 0, 0, 1, t10; 0, 0, 0, 1; t8 * t20, -t8 * t18, t6, t8 * pkin(2) + t6 * qJ(3) + t25; t6 * t20, -t6 * t18, -t8, t6 * pkin(2) - t8 * qJ(3) + t26; t18, t20, 0, t10; 0, 0, 0, 1; t8 * t7, -t8 * t5, t6, -t6 * t22 + t8 * t4 + t25; t6 * t7, -t6 * t5, -t8, t8 * t22 + t6 * t4 + t26; t5, t7, 0, t23; 0, 0, 0, 1; t8 * t3, -t8 * t2, t6, t8 * t1 - t6 * t15 + t25; t6 * t3, -t6 * t2, -t8, t6 * t1 + t8 * t15 + t26; t2, t3, 0, pkin(4) * t5 + t23; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
