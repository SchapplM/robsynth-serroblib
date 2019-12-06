% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:02:46
% EndTime: 2019-12-05 17:02:46
% DurationCPUTime: 0.08s
% Computational Cost: add. (46->22), mult. (34->21), div. (0->0), fcn. (67->8), ass. (0->23)
t11 = sin(qJ(2));
t8 = qJ(3) + qJ(4);
t4 = sin(t8);
t22 = t11 * t4;
t9 = sin(qJ(5));
t21 = t11 * t9;
t14 = cos(qJ(2));
t20 = t14 * t4;
t19 = t14 * t9;
t12 = cos(qJ(5));
t18 = t11 * t12;
t13 = cos(qJ(3));
t17 = t11 * t13;
t16 = t14 * t12;
t15 = t14 * t13;
t7 = pkin(1) + 0;
t6 = qJ(1) + 0;
t10 = sin(qJ(3));
t5 = cos(t8);
t3 = -t10 * pkin(2) + 0;
t2 = pkin(2) * t15 + t7;
t1 = pkin(2) * t17 + t6;
t23 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, t6; 0, 0, 0, 1; t14, -t11, 0, t7; 0, 0, -1, 0; t11, t14, 0, t6; 0, 0, 0, 1; t15, -t14 * t10, t11, t7; -t10, -t13, 0, 0; t17, -t11 * t10, -t14, t6; 0, 0, 0, 1; t14 * t5, -t20, t11, t2; -t4, -t5, 0, t3; t11 * t5, -t22, -t14, t1; 0, 0, 0, 1; t5 * t16 + t21, -t5 * t19 + t18, t20, t2; -t4 * t12, t4 * t9, t5, t3; t5 * t18 - t19, -t5 * t21 - t16, t22, t1; 0, 0, 0, 1;];
T_ges = t23;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
