% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-10-24 10:51
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRRP2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:51:37
% EndTime: 2019-10-24 10:51:37
% DurationCPUTime: 0.08s
% Computational Cost: add. (105->36), mult. (41->26), div. (0->0), fcn. (73->8), ass. (0->24)
t15 = qJ(1) + qJ(2);
t6 = sin(t15);
t14 = qJ(3) + qJ(4);
t7 = cos(t14);
t26 = t6 * t7;
t5 = sin(t14);
t8 = cos(t15);
t25 = t8 * t5;
t20 = -pkin(8) - pkin(7);
t18 = cos(qJ(3));
t4 = t18 * pkin(3) + pkin(2);
t24 = pkin(5) + 0;
t19 = cos(qJ(1));
t23 = t19 * pkin(1) + 0;
t9 = pkin(6) + t24;
t17 = sin(qJ(1));
t22 = -t17 * pkin(1) + 0;
t16 = sin(qJ(3));
t21 = t16 * pkin(3) + t9;
t13 = -qJ(5) + t20;
t3 = t8 * t7;
t2 = t6 * t5;
t1 = pkin(4) * t7 + t4;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t24; -t17, -t19, 0, 0; t19, -t17, 0, 0; 0, 0, 0, 1; 0, 0, 1, t9; -t6, -t8, 0, t22; t8, -t6, 0, t23; 0, 0, 0, 1; t16, t18, 0, t9; -t6 * t18, t6 * t16, t8, -t6 * pkin(2) + t8 * pkin(7) + t22; t8 * t18, -t8 * t16, t6, t8 * pkin(2) + t6 * pkin(7) + t23; 0, 0, 0, 1; t5, t7, 0, t21; -t26, t2, t8, -t8 * t20 - t6 * t4 + t22; t3, -t25, t6, -t6 * t20 + t8 * t4 + t23; 0, 0, 0, 1; t5, t7, 0, pkin(4) * t5 + t21; -t26, t2, t8, -t6 * t1 - t8 * t13 + t22; t3, -t25, t6, t8 * t1 - t6 * t13 + t23; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
