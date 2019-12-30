% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4PRRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 12:18:23
% EndTime: 2019-12-29 12:18:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (52->32), mult. (75->36), div. (0->0), fcn. (113->6), ass. (0->26)
t13 = sin(pkin(6));
t16 = sin(qJ(3));
t28 = t13 * t16;
t19 = cos(qJ(2));
t27 = t16 * t19;
t17 = sin(qJ(2));
t26 = t17 * t16;
t18 = cos(qJ(3));
t25 = t18 * t19;
t24 = t13 * pkin(1) + 0;
t12 = qJ(1) + 0;
t14 = cos(pkin(6));
t23 = t14 * pkin(1) + t13 * pkin(4) + 0;
t22 = pkin(2) * t19 + pkin(5) * t17;
t15 = -qJ(4) - pkin(5);
t8 = t18 * pkin(3) + pkin(2);
t21 = -t15 * t17 + t19 * t8;
t20 = -t14 * pkin(4) + t24;
t7 = t17 * t18;
t6 = t14 * t17;
t5 = t13 * t17;
t4 = t14 * t25 + t28;
t3 = t13 * t18 - t14 * t27;
t2 = t13 * t25 - t14 * t16;
t1 = -t13 * t27 - t14 * t18;
t9 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t14, -t13, 0, 0; t13, t14, 0, 0; 0, 0, 1, t12; 0, 0, 0, 1; t14 * t19, -t6, t13, t23; t13 * t19, -t5, -t14, t20; t17, t19, 0, t12; 0, 0, 0, 1; t4, t3, t6, t22 * t14 + t23; t2, t1, t5, t22 * t13 + t20; t7, -t26, -t19, t17 * pkin(2) - t19 * pkin(5) + t12; 0, 0, 0, 1; t4, t3, t6, pkin(3) * t28 + t21 * t14 + t23; t2, t1, t5, (-pkin(3) * t16 - pkin(4)) * t14 + t21 * t13 + t24; t7, -t26, -t19, t19 * t15 + t17 * t8 + t12; 0, 0, 0, 1;];
T_ges = t9;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
