% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR9
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
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRPR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:05
% EndTime: 2019-12-31 21:22:06
% DurationCPUTime: 0.14s
% Computational Cost: add. (120->55), mult. (117->64), div. (0->0), fcn. (168->10), ass. (0->31)
t19 = sin(qJ(3));
t34 = t19 * pkin(3);
t22 = cos(qJ(3));
t7 = t22 * pkin(3) + pkin(2);
t21 = sin(qJ(1));
t33 = t21 * t19;
t23 = cos(qJ(2));
t32 = t21 * t23;
t24 = cos(qJ(1));
t31 = t24 * t23;
t18 = -qJ(4) - pkin(7);
t17 = pkin(5) + 0;
t16 = qJ(3) + pkin(9);
t30 = t21 * pkin(1) + 0;
t29 = t24 * pkin(1) + t21 * pkin(6) + 0;
t20 = sin(qJ(2));
t28 = pkin(2) * t23 + pkin(7) * t20;
t9 = cos(t16);
t1 = pkin(4) * t9 + t7;
t15 = -pkin(8) + t18;
t27 = t1 * t23 - t15 * t20;
t26 = -t18 * t20 + t23 * t7;
t25 = -t24 * pkin(6) + t30;
t10 = qJ(5) + t16;
t8 = sin(t16);
t6 = t24 * t20;
t5 = t21 * t20;
t4 = cos(t10);
t3 = sin(t10);
t2 = pkin(4) * t8 + t34;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t21, 0, 0; t21, t24, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; t31, -t6, t21, t29; t32, -t5, -t24, t25; t20, t23, 0, t17; 0, 0, 0, 1; t22 * t31 + t33, -t19 * t31 + t21 * t22, t6, t28 * t24 + t29; -t24 * t19 + t22 * t32, -t19 * t32 - t24 * t22, t5, t28 * t21 + t25; t20 * t22, -t20 * t19, -t23, t20 * pkin(2) - t23 * pkin(7) + t17; 0, 0, 0, 1; t21 * t8 + t9 * t31, t21 * t9 - t8 * t31, t6, pkin(3) * t33 + t26 * t24 + t29; -t24 * t8 + t9 * t32, -t24 * t9 - t8 * t32, t5, (-pkin(6) - t34) * t24 + t26 * t21 + t30; t20 * t9, -t20 * t8, -t23, t23 * t18 + t20 * t7 + t17; 0, 0, 0, 1; t21 * t3 + t4 * t31, t21 * t4 - t3 * t31, t6, t21 * t2 + t27 * t24 + t29; -t24 * t3 + t4 * t32, -t24 * t4 - t3 * t32, t5, (-pkin(6) - t2) * t24 + t27 * t21 + t30; t20 * t4, -t20 * t3, -t23, t20 * t1 + t23 * t15 + t17; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
