% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-10-24 10:46
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:46:12
% EndTime: 2019-10-24 10:46:12
% DurationCPUTime: 0.11s
% Computational Cost: add. (111->37), mult. (41->30), div. (0->0), fcn. (73->10), ass. (0->23)
t21 = -pkin(7) - pkin(6);
t19 = cos(qJ(3));
t4 = t19 * pkin(3) + pkin(2);
t16 = qJ(3) + qJ(4);
t25 = pkin(5) + 0;
t20 = cos(qJ(1));
t24 = t20 * pkin(1) + 0;
t7 = qJ(2) + t25;
t18 = sin(qJ(1));
t23 = -t18 * pkin(1) + 0;
t17 = sin(qJ(3));
t22 = t17 * pkin(3) + t7;
t15 = -pkin(8) + t21;
t14 = qJ(1) + pkin(9);
t10 = qJ(5) + t16;
t9 = cos(t16);
t8 = sin(t16);
t6 = cos(t14);
t5 = sin(t14);
t3 = cos(t10);
t2 = sin(t10);
t1 = pkin(4) * t9 + t4;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t25; -t18, -t20, 0, 0; t20, -t18, 0, 0; 0, 0, 0, 1; 0, 0, 1, t7; -t5, -t6, 0, t23; t6, -t5, 0, t24; 0, 0, 0, 1; t17, t19, 0, t7; -t5 * t19, t5 * t17, t6, -t5 * pkin(2) + t6 * pkin(6) + t23; t6 * t19, -t6 * t17, t5, t6 * pkin(2) + t5 * pkin(6) + t24; 0, 0, 0, 1; t8, t9, 0, t22; -t5 * t9, t5 * t8, t6, -t6 * t21 - t5 * t4 + t23; t6 * t9, -t6 * t8, t5, -t5 * t21 + t6 * t4 + t24; 0, 0, 0, 1; t2, t3, 0, pkin(4) * t8 + t22; -t5 * t3, t5 * t2, t6, -t5 * t1 - t6 * t15 + t23; t6 * t3, -t6 * t2, t5, t6 * t1 - t5 * t15 + t24; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
