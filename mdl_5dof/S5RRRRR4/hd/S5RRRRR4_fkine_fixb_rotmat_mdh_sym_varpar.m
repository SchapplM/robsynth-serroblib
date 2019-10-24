% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-10-24 10:52
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:52:00
% EndTime: 2019-10-24 10:52:00
% DurationCPUTime: 0.09s
% Computational Cost: add. (112->42), mult. (50->38), div. (0->0), fcn. (86->10), ass. (0->23)
t24 = -pkin(7) - pkin(6);
t22 = cos(qJ(2));
t9 = t22 * pkin(2) + pkin(1);
t19 = qJ(2) + qJ(3);
t18 = -pkin(8) + t24;
t17 = pkin(5) + 0;
t12 = cos(t19);
t2 = pkin(3) * t12 + t9;
t13 = qJ(4) + t19;
t20 = sin(qJ(2));
t26 = t20 * pkin(2) + t17;
t11 = sin(t19);
t25 = pkin(3) * t11 + t26;
t23 = cos(qJ(1));
t21 = sin(qJ(1));
t14 = -pkin(9) + t18;
t10 = qJ(5) + t13;
t8 = cos(t13);
t7 = sin(t13);
t4 = cos(t10);
t3 = sin(t10);
t1 = pkin(4) * t8 + t2;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t21, 0, 0; t21, t23, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; t23 * t22, -t23 * t20, t21, t23 * pkin(1) + t21 * pkin(6) + 0; t21 * t22, -t21 * t20, -t23, t21 * pkin(1) - t23 * pkin(6) + 0; t20, t22, 0, t17; 0, 0, 0, 1; t23 * t12, -t23 * t11, t21, -t21 * t24 + t23 * t9 + 0; t21 * t12, -t21 * t11, -t23, t21 * t9 + t23 * t24 + 0; t11, t12, 0, t26; 0, 0, 0, 1; t23 * t8, -t23 * t7, t21, -t21 * t18 + t23 * t2 + 0; t21 * t8, -t21 * t7, -t23, t23 * t18 + t21 * t2 + 0; t7, t8, 0, t25; 0, 0, 0, 1; t23 * t4, -t23 * t3, t21, t23 * t1 - t21 * t14 + 0; t21 * t4, -t21 * t3, -t23, t21 * t1 + t23 * t14 + 0; t3, t4, 0, pkin(4) * t7 + t25; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
