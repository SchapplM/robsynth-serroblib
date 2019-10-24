% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR6
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

function T_c_mdh = S5RRRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:52:46
% EndTime: 2019-10-24 10:52:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (111->37), mult. (41->30), div. (0->0), fcn. (73->10), ass. (0->23)
t21 = -pkin(8) - pkin(7);
t19 = cos(qJ(3));
t4 = t19 * pkin(3) + pkin(2);
t15 = qJ(3) + qJ(4);
t25 = pkin(5) + 0;
t20 = cos(qJ(1));
t24 = t20 * pkin(1) + 0;
t10 = pkin(6) + t25;
t18 = sin(qJ(1));
t23 = -t18 * pkin(1) + 0;
t17 = sin(qJ(3));
t22 = t17 * pkin(3) + t10;
t16 = qJ(1) + qJ(2);
t14 = -pkin(9) + t21;
t9 = qJ(5) + t15;
t8 = cos(t16);
t7 = cos(t15);
t6 = sin(t16);
t5 = sin(t15);
t3 = cos(t9);
t2 = sin(t9);
t1 = pkin(4) * t7 + t4;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t25; -t18, -t20, 0, 0; t20, -t18, 0, 0; 0, 0, 0, 1; 0, 0, 1, t10; -t6, -t8, 0, t23; t8, -t6, 0, t24; 0, 0, 0, 1; t17, t19, 0, t10; -t6 * t19, t6 * t17, t8, -t6 * pkin(2) + t8 * pkin(7) + t23; t8 * t19, -t8 * t17, t6, t8 * pkin(2) + t6 * pkin(7) + t24; 0, 0, 0, 1; t5, t7, 0, t22; -t6 * t7, t6 * t5, t8, -t8 * t21 - t6 * t4 + t23; t8 * t7, -t8 * t5, t6, -t6 * t21 + t8 * t4 + t24; 0, 0, 0, 1; t2, t3, 0, pkin(4) * t5 + t22; -t6 * t3, t6 * t2, t8, -t6 * t1 - t8 * t14 + t23; t8 * t3, -t8 * t2, t6, t8 * t1 - t6 * t14 + t24; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
