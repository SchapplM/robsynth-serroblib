% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-29 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRPR14_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 17:06:51
% EndTime: 2019-12-29 17:06:51
% DurationCPUTime: 0.24s
% Computational Cost: add. (83->38), mult. (66->36), div. (0->0), fcn. (103->8), ass. (0->28)
t13 = sin(qJ(1));
t8 = qJ(3) + pkin(8);
t3 = cos(t8);
t32 = t13 * t3;
t16 = cos(qJ(1));
t31 = t16 * t3;
t11 = sin(qJ(5));
t30 = t13 * t11;
t12 = sin(qJ(3));
t29 = t13 * t12;
t14 = cos(qJ(5));
t28 = t13 * t14;
t27 = t16 * t11;
t26 = t16 * t14;
t9 = pkin(5) + 0;
t25 = t13 * pkin(1) + 0;
t24 = pkin(2) + t9;
t23 = t16 * pkin(1) + t13 * qJ(2) + 0;
t22 = -pkin(3) * t12 - qJ(2);
t15 = cos(qJ(3));
t21 = t15 * pkin(3) + t24;
t2 = sin(t8);
t20 = pkin(4) * t2 - pkin(7) * t3;
t10 = -qJ(4) - pkin(6);
t19 = -t13 * t10 + t25;
t18 = -t16 * qJ(2) + t25;
t17 = pkin(3) * t29 - t16 * t10 + t23;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t13, 0, 0; t13, t16, 0, 0; 0, 0, 1, t9; 0, 0, 0, 1; 0, -t16, t13, t23; 0, -t13, -t16, t18; 1, 0, 0, t9; 0, 0, 0, 1; t29, t13 * t15, t16, t16 * pkin(6) + t23; -t16 * t12, -t16 * t15, t13, t13 * pkin(6) + t18; t15, -t12, 0, t24; 0, 0, 0, 1; t13 * t2, t32, t16, t17; -t16 * t2, -t31, t13, t22 * t16 + t19; t3, -t2, 0, t21; 0, 0, 0, 1; t2 * t28 + t27, -t2 * t30 + t26, -t32, t20 * t13 + t17; -t2 * t26 + t30, t2 * t27 + t28, t31, (-t20 + t22) * t16 + t19; t3 * t14, -t3 * t11, t2, t3 * pkin(4) + t2 * pkin(7) + t21; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
