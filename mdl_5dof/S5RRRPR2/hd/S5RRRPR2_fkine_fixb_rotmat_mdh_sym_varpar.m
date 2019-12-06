% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:39
% EndTime: 2019-12-05 18:40:39
% DurationCPUTime: 0.08s
% Computational Cost: add. (118->31), mult. (26->14), div. (0->0), fcn. (50->10), ass. (0->24)
t13 = qJ(1) + qJ(2);
t26 = pkin(5) + 0;
t17 = cos(qJ(1));
t25 = t17 * pkin(1) + 0;
t24 = pkin(6) + t26;
t11 = qJ(3) + t13;
t10 = cos(t13);
t23 = pkin(2) * t10 + t25;
t22 = pkin(7) + t24;
t15 = sin(qJ(1));
t21 = -t15 * pkin(1) + 0;
t8 = cos(t11);
t20 = pkin(3) * t8 + t23;
t9 = sin(t13);
t19 = -pkin(2) * t9 + t21;
t7 = sin(t11);
t18 = -pkin(3) * t7 + t19;
t16 = cos(qJ(5));
t14 = sin(qJ(5));
t6 = pkin(9) + t11;
t4 = qJ(4) + t22;
t2 = cos(t6);
t1 = sin(t6);
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t26; -t15, -t17, 0, 0; t17, -t15, 0, 0; 0, 0, 0, 1; 0, 0, 1, t24; -t9, -t10, 0, t21; t10, -t9, 0, t25; 0, 0, 0, 1; 0, 0, 1, t22; -t7, -t8, 0, t19; t8, -t7, 0, t23; 0, 0, 0, 1; 0, 0, 1, t4; -t1, -t2, 0, t18; t2, -t1, 0, t20; 0, 0, 0, 1; t14, t16, 0, t4; -t1 * t16, t1 * t14, t2, -t1 * pkin(4) + t2 * pkin(8) + t18; t2 * t16, -t2 * t14, t1, t2 * pkin(4) + t1 * pkin(8) + t20; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
