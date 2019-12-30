% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-29 20:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 20:01:06
% EndTime: 2019-12-29 20:01:07
% DurationCPUTime: 0.21s
% Computational Cost: add. (108->39), mult. (84->38), div. (0->0), fcn. (126->8), ass. (0->27)
t22 = cos(qJ(1));
t16 = qJ(2) + qJ(3);
t12 = sin(t16);
t32 = qJ(4) * t12;
t13 = cos(t16);
t8 = t22 * t13;
t35 = pkin(3) * t8 + t22 * t32;
t19 = sin(qJ(1));
t34 = t19 * t12;
t7 = t19 * t13;
t33 = t22 * t12;
t15 = pkin(5) + 0;
t21 = cos(qJ(2));
t10 = t21 * pkin(2) + pkin(1);
t31 = t22 * t10 + 0;
t18 = sin(qJ(2));
t30 = t18 * pkin(2) + t15;
t23 = -pkin(7) - pkin(6);
t29 = t19 * t10 + t22 * t23 + 0;
t28 = pkin(3) * t7 + t19 * t32 + t29;
t17 = sin(qJ(5));
t20 = cos(qJ(5));
t27 = t12 * t20 - t13 * t17;
t26 = t12 * t17 + t13 * t20;
t25 = -t19 * t23 + t31;
t24 = t12 * pkin(3) - t13 * qJ(4) + t30;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t19, 0, 0; t19, t22, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t22 * t21, -t22 * t18, t19, t22 * pkin(1) + t19 * pkin(6) + 0; t19 * t21, -t19 * t18, -t22, t19 * pkin(1) - t22 * pkin(6) + 0; t18, t21, 0, t15; 0, 0, 0, 1; t8, -t33, t19, t25; t7, -t34, -t22, t29; t12, t13, 0, t30; 0, 0, 0, 1; t8, t19, t33, t25 + t35; t7, -t22, t34, t28; t12, 0, -t13, t24; 0, 0, 0, 1; t26 * t22, t27 * t22, -t19, pkin(4) * t8 + (-pkin(8) - t23) * t19 + t31 + t35; t26 * t19, t27 * t19, t22, pkin(4) * t7 + t22 * pkin(8) + t28; t27, -t26, 0, t12 * pkin(4) + t24; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
