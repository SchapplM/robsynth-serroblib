% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:31
% EndTime: 2019-12-05 18:08:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (31->21), mult. (65->28), div. (0->0), fcn. (115->8), ass. (0->21)
t8 = sin(qJ(4));
t9 = sin(qJ(3));
t20 = t9 * t8;
t10 = sin(qJ(1));
t19 = t10 * t9;
t14 = cos(qJ(1));
t18 = t14 * t9;
t12 = cos(qJ(4));
t17 = t9 * t12;
t13 = cos(qJ(3));
t16 = t10 * t13;
t15 = t14 * t13;
t11 = cos(qJ(5));
t7 = sin(qJ(5));
t6 = -t14 * qJ(2) + 0;
t5 = t10 * qJ(2) + 0;
t4 = t10 * t8 + t12 * t15;
t3 = -t10 * t12 + t8 * t15;
t2 = t12 * t16 - t14 * t8;
t1 = t14 * t12 + t8 * t16;
t21 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t14, -t10, 0, 0; t10, t14, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t14, 0, t10, t5; t10, 0, -t14, t6; 0, 1, 0, 0; 0, 0, 0, 1; t15, -t18, t10, t5; t16, -t19, -t14, t6; t9, t13, 0, 0; 0, 0, 0, 1; t4, -t3, t18, t5; t2, -t1, t19, t6; t17, -t20, -t13, 0; 0, 0, 0, 1; t4 * t11 + t7 * t18, t11 * t18 - t4 * t7, t3, t5; t2 * t11 + t7 * t19, t11 * t19 - t2 * t7, t1, t6; t11 * t17 - t13 * t7, -t13 * t11 - t7 * t17, t20, 0; 0, 0, 0, 1;];
T_ges = t21;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
