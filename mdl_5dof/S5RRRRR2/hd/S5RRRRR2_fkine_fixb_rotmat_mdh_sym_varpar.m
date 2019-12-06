% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR2
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:50
% EndTime: 2019-12-05 18:52:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (73->21), mult. (42->23), div. (0->0), fcn. (79->10), ass. (0->26)
t12 = qJ(3) + qJ(4);
t6 = sin(t12);
t13 = qJ(1) + qJ(2);
t7 = sin(t13);
t27 = t7 * t6;
t9 = cos(t13);
t26 = t9 * t6;
t14 = sin(qJ(5));
t25 = t7 * t14;
t17 = cos(qJ(5));
t24 = t7 * t17;
t18 = cos(qJ(3));
t23 = t7 * t18;
t22 = t9 * t14;
t21 = t9 * t17;
t20 = t9 * t18;
t16 = sin(qJ(1));
t4 = t16 * pkin(1) + 0;
t19 = cos(qJ(1));
t5 = t19 * pkin(1) + 0;
t15 = sin(qJ(3));
t8 = cos(t12);
t3 = t15 * pkin(2) + 0;
t2 = pkin(2) * t20 + t5;
t1 = pkin(2) * t23 + t4;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t16, 0, 0; t16, t19, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t9, -t7, 0, t5; t7, t9, 0, t4; 0, 0, 1, 0; 0, 0, 0, 1; t20, -t9 * t15, t7, t5; t23, -t7 * t15, -t9, t4; t15, t18, 0, 0; 0, 0, 0, 1; t9 * t8, -t26, t7, t2; t7 * t8, -t27, -t9, t1; t6, t8, 0, t3; 0, 0, 0, 1; t8 * t21 + t25, -t8 * t22 + t24, t26, t2; t8 * t24 - t22, -t8 * t25 - t21, t27, t1; t6 * t17, -t6 * t14, -t8, t3; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
