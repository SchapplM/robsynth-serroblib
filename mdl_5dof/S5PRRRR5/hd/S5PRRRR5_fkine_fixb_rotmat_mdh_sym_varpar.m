% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-10-24 10:36
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:36:55
% EndTime: 2019-10-24 10:36:55
% DurationCPUTime: 0.09s
% Computational Cost: add. (110->36), mult. (41->30), div. (0->0), fcn. (73->10), ass. (0->23)
t22 = -pkin(7) - pkin(6);
t21 = cos(qJ(3));
t4 = t21 * pkin(3) + pkin(2);
t17 = qJ(3) + qJ(4);
t18 = sin(pkin(9));
t26 = t18 * pkin(1) + 0;
t19 = cos(pkin(9));
t25 = t19 * pkin(1) + 0;
t24 = qJ(1) + 0;
t7 = pkin(5) + t24;
t20 = sin(qJ(3));
t23 = t20 * pkin(3) + t7;
t16 = -pkin(8) + t22;
t15 = pkin(9) + qJ(2);
t12 = qJ(5) + t17;
t9 = cos(t17);
t8 = sin(t17);
t6 = cos(t15);
t5 = sin(t15);
t3 = cos(t12);
t2 = sin(t12);
t1 = pkin(4) * t9 + t4;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t18, 0, 0; t18, t19, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t6, -t5, 0, t25; t5, t6, 0, t26; 0, 0, 1, t7; 0, 0, 0, 1; t6 * t21, -t6 * t20, t5, t6 * pkin(2) + t5 * pkin(6) + t25; t5 * t21, -t5 * t20, -t6, t5 * pkin(2) - t6 * pkin(6) + t26; t20, t21, 0, t7; 0, 0, 0, 1; t6 * t9, -t6 * t8, t5, -t5 * t22 + t6 * t4 + t25; t5 * t9, -t5 * t8, -t6, t6 * t22 + t5 * t4 + t26; t8, t9, 0, t23; 0, 0, 0, 1; t6 * t3, -t6 * t2, t5, t6 * t1 - t5 * t16 + t25; t5 * t3, -t5 * t2, -t6, t5 * t1 + t6 * t16 + t26; t2, t3, 0, pkin(4) * t8 + t23; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
