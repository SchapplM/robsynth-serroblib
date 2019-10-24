% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP1
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

function T_c_mdh = S5RRRRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:51:16
% EndTime: 2019-10-24 10:51:16
% DurationCPUTime: 0.09s
% Computational Cost: add. (106->41), mult. (50->34), div. (0->0), fcn. (86->8), ass. (0->24)
t23 = -pkin(7) - pkin(6);
t20 = sin(qJ(1));
t18 = qJ(2) + qJ(3);
t13 = qJ(4) + t18;
t7 = sin(t13);
t27 = t20 * t7;
t22 = cos(qJ(1));
t26 = t22 * t7;
t21 = cos(qJ(2));
t9 = t21 * pkin(2) + pkin(1);
t17 = -pkin(8) + t23;
t16 = pkin(5) + 0;
t11 = cos(t18);
t2 = pkin(3) * t11 + t9;
t19 = sin(qJ(2));
t25 = t19 * pkin(2) + t16;
t10 = sin(t18);
t24 = pkin(3) * t10 + t25;
t12 = -qJ(5) + t17;
t8 = cos(t13);
t4 = t22 * t8;
t3 = t20 * t8;
t1 = pkin(4) * t8 + t2;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t20, 0, 0; t20, t22, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t22 * t21, -t22 * t19, t20, t22 * pkin(1) + t20 * pkin(6) + 0; t20 * t21, -t20 * t19, -t22, t20 * pkin(1) - t22 * pkin(6) + 0; t19, t21, 0, t16; 0, 0, 0, 1; t22 * t11, -t22 * t10, t20, -t20 * t23 + t22 * t9 + 0; t20 * t11, -t20 * t10, -t22, t20 * t9 + t22 * t23 + 0; t10, t11, 0, t25; 0, 0, 0, 1; t4, -t26, t20, -t20 * t17 + t22 * t2 + 0; t3, -t27, -t22, t22 * t17 + t20 * t2 + 0; t7, t8, 0, t24; 0, 0, 0, 1; t4, -t26, t20, t22 * t1 - t20 * t12 + 0; t3, -t27, -t22, t20 * t1 + t22 * t12 + 0; t7, t8, 0, pkin(4) * t7 + t24; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
