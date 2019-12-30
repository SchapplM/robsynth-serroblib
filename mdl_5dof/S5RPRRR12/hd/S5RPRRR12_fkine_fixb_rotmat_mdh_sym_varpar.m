% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-29 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 17:58:50
% EndTime: 2019-12-29 17:58:50
% DurationCPUTime: 0.19s
% Computational Cost: add. (83->38), mult. (66->36), div. (0->0), fcn. (103->8), ass. (0->28)
t12 = sin(qJ(1));
t9 = qJ(3) + qJ(4);
t3 = cos(t9);
t32 = t12 * t3;
t15 = cos(qJ(1));
t31 = t15 * t3;
t10 = sin(qJ(5));
t30 = t12 * t10;
t11 = sin(qJ(3));
t29 = t12 * t11;
t13 = cos(qJ(5));
t28 = t12 * t13;
t27 = t15 * t10;
t26 = t15 * t13;
t8 = pkin(5) + 0;
t25 = t12 * pkin(1) + 0;
t24 = pkin(2) + t8;
t23 = t15 * pkin(1) + t12 * qJ(2) + 0;
t22 = -pkin(3) * t11 - qJ(2);
t14 = cos(qJ(3));
t21 = t14 * pkin(3) + t24;
t2 = sin(t9);
t20 = pkin(4) * t2 - pkin(8) * t3;
t16 = -pkin(7) - pkin(6);
t19 = -t12 * t16 + t25;
t18 = -t15 * qJ(2) + t25;
t17 = pkin(3) * t29 - t15 * t16 + t23;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t15, -t12, 0, 0; t12, t15, 0, 0; 0, 0, 1, t8; 0, 0, 0, 1; 0, -t15, t12, t23; 0, -t12, -t15, t18; 1, 0, 0, t8; 0, 0, 0, 1; t29, t12 * t14, t15, t15 * pkin(6) + t23; -t15 * t11, -t15 * t14, t12, t12 * pkin(6) + t18; t14, -t11, 0, t24; 0, 0, 0, 1; t12 * t2, t32, t15, t17; -t15 * t2, -t31, t12, t22 * t15 + t19; t3, -t2, 0, t21; 0, 0, 0, 1; t2 * t28 + t27, -t2 * t30 + t26, -t32, t20 * t12 + t17; -t2 * t26 + t30, t2 * t27 + t28, t31, (-t20 + t22) * t15 + t19; t3 * t13, -t3 * t10, t2, t3 * pkin(4) + t2 * pkin(8) + t21; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
