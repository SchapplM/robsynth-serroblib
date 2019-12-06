% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR2
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:26:57
% EndTime: 2019-12-05 18:26:57
% DurationCPUTime: 0.09s
% Computational Cost: add. (112->42), mult. (50->38), div. (0->0), fcn. (86->10), ass. (0->23)
t23 = cos(qJ(2));
t10 = t23 * pkin(2) + pkin(1);
t20 = -qJ(3) - pkin(6);
t19 = pkin(5) + 0;
t18 = qJ(2) + pkin(9);
t12 = cos(t18);
t2 = pkin(3) * t12 + t10;
t17 = -pkin(7) + t20;
t21 = sin(qJ(2));
t26 = t21 * pkin(2) + t19;
t13 = qJ(4) + t18;
t11 = sin(t18);
t25 = pkin(3) * t11 + t26;
t24 = cos(qJ(1));
t22 = sin(qJ(1));
t14 = -pkin(8) + t17;
t9 = qJ(5) + t13;
t8 = cos(t13);
t7 = sin(t13);
t4 = cos(t9);
t3 = sin(t9);
t1 = pkin(4) * t8 + t2;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t22, 0, 0; t22, t24, 0, 0; 0, 0, 1, t19; 0, 0, 0, 1; t24 * t23, -t24 * t21, t22, t24 * pkin(1) + t22 * pkin(6) + 0; t22 * t23, -t22 * t21, -t24, t22 * pkin(1) - t24 * pkin(6) + 0; t21, t23, 0, t19; 0, 0, 0, 1; t24 * t12, -t24 * t11, t22, t24 * t10 - t22 * t20 + 0; t22 * t12, -t22 * t11, -t24, t22 * t10 + t24 * t20 + 0; t11, t12, 0, t26; 0, 0, 0, 1; t24 * t8, -t24 * t7, t22, -t22 * t17 + t24 * t2 + 0; t22 * t8, -t22 * t7, -t24, t24 * t17 + t22 * t2 + 0; t7, t8, 0, t25; 0, 0, 0, 1; t24 * t4, -t24 * t3, t22, t24 * t1 - t22 * t14 + 0; t22 * t4, -t22 * t3, -t24, t22 * t1 + t24 * t14 + 0; t3, t4, 0, pkin(4) * t7 + t25; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
