% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:18:56
% EndTime: 2019-12-05 16:18:56
% DurationCPUTime: 0.09s
% Computational Cost: add. (110->36), mult. (41->30), div. (0->0), fcn. (73->10), ass. (0->23)
t22 = cos(qJ(3));
t4 = t22 * pkin(3) + pkin(2);
t20 = -qJ(4) - pkin(6);
t17 = qJ(3) + pkin(9);
t18 = sin(pkin(8));
t26 = t18 * pkin(1) + 0;
t19 = cos(pkin(8));
t25 = t19 * pkin(1) + 0;
t24 = qJ(1) + 0;
t10 = pkin(5) + t24;
t21 = sin(qJ(3));
t23 = t21 * pkin(3) + t10;
t16 = pkin(8) + qJ(2);
t15 = -pkin(7) + t20;
t9 = qJ(5) + t17;
t8 = cos(t17);
t7 = cos(t16);
t6 = sin(t17);
t5 = sin(t16);
t3 = cos(t9);
t2 = sin(t9);
t1 = pkin(4) * t8 + t4;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t18, 0, 0; t18, t19, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t7, -t5, 0, t25; t5, t7, 0, t26; 0, 0, 1, t10; 0, 0, 0, 1; t7 * t22, -t7 * t21, t5, t7 * pkin(2) + t5 * pkin(6) + t25; t5 * t22, -t5 * t21, -t7, t5 * pkin(2) - t7 * pkin(6) + t26; t21, t22, 0, t10; 0, 0, 0, 1; t7 * t8, -t7 * t6, t5, -t5 * t20 + t7 * t4 + t25; t5 * t8, -t5 * t6, -t7, t7 * t20 + t5 * t4 + t26; t6, t8, 0, t23; 0, 0, 0, 1; t7 * t3, -t7 * t2, t5, t7 * t1 - t5 * t15 + t25; t5 * t3, -t5 * t2, -t7, t5 * t1 + t7 * t15 + t26; t2, t3, 0, pkin(4) * t6 + t23; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
