% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR3
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
% Datum: 2019-10-24 10:45
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:45:08
% EndTime: 2019-10-24 10:45:08
% DurationCPUTime: 0.09s
% Computational Cost: add. (79->38), mult. (47->30), div. (0->0), fcn. (79->8), ass. (0->22)
t18 = -pkin(7) - pkin(6);
t14 = sin(qJ(3));
t25 = t14 * pkin(3);
t15 = sin(qJ(1));
t24 = t15 * t14;
t13 = qJ(3) + qJ(4);
t11 = pkin(5) + 0;
t23 = t15 * pkin(1) + 0;
t22 = pkin(2) + t11;
t17 = cos(qJ(1));
t21 = t17 * pkin(1) + t15 * qJ(2) + 0;
t16 = cos(qJ(3));
t20 = t16 * pkin(3) + t22;
t19 = -t17 * qJ(2) + t23;
t12 = -pkin(8) + t18;
t7 = qJ(5) + t13;
t5 = cos(t13);
t4 = sin(t13);
t3 = cos(t7);
t2 = sin(t7);
t1 = pkin(4) * t4 + t25;
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t15, 0, 0; t15, t17, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; 0, -t17, t15, t21; 0, -t15, -t17, t19; 1, 0, 0, t11; 0, 0, 0, 1; t24, t15 * t16, t17, t17 * pkin(6) + t21; -t17 * t14, -t17 * t16, t15, t15 * pkin(6) + t19; t16, -t14, 0, t22; 0, 0, 0, 1; t15 * t4, t15 * t5, t17, pkin(3) * t24 - t17 * t18 + t21; -t17 * t4, -t17 * t5, t15, -t15 * t18 + (-qJ(2) - t25) * t17 + t23; t5, -t4, 0, t20; 0, 0, 0, 1; t15 * t2, t15 * t3, t17, t15 * t1 - t17 * t12 + t21; -t17 * t2, -t17 * t3, t15, -t15 * t12 + (-qJ(2) - t1) * t17 + t23; t3, -t2, 0, pkin(4) * t5 + t20; 0, 0, 0, 1;];
T_ges = t6;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
