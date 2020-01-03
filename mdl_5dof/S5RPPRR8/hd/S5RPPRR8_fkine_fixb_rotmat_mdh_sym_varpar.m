% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:00:59
% EndTime: 2019-12-31 18:00:59
% DurationCPUTime: 0.09s
% Computational Cost: add. (96->35), mult. (78->28), div. (0->0), fcn. (120->8), ass. (0->24)
t17 = sin(pkin(8));
t20 = sin(qJ(1));
t32 = t20 * t17;
t16 = pkin(5) + 0;
t31 = pkin(8) + qJ(4);
t30 = t20 * pkin(1) + 0;
t29 = -qJ(3) + t16;
t22 = cos(qJ(1));
t28 = t22 * pkin(1) + t20 * qJ(2) + 0;
t27 = cos(t31);
t26 = sin(t31);
t18 = cos(pkin(8));
t11 = t18 * pkin(3) + pkin(2);
t25 = pkin(3) * t32 + t22 * t11 + t28;
t24 = -t22 * qJ(2) + t30;
t23 = t20 * t11 + (-pkin(3) * t17 - qJ(2)) * t22 + t30;
t21 = cos(qJ(5));
t19 = sin(qJ(5));
t12 = -pkin(6) + t29;
t4 = -t22 * t17 + t20 * t18;
t3 = -t22 * t18 - t32;
t2 = -t20 * t27 + t22 * t26;
t1 = -t20 * t26 - t22 * t27;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t20, 0, 0; t20, t22, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t22, 0, t20, t28; t20, 0, -t22, t24; 0, 1, 0, t16; 0, 0, 0, 1; -t3, t4, 0, t22 * pkin(2) + t28; t4, t3, 0, t20 * pkin(2) + t24; 0, 0, -1, t29; 0, 0, 0, 1; -t1, -t2, 0, t25; -t2, t1, 0, t23; 0, 0, -1, t12; 0, 0, 0, 1; -t1 * t21, t1 * t19, t2, -t1 * pkin(4) + t2 * pkin(7) + t25; -t2 * t21, t2 * t19, -t1, -t2 * pkin(4) - t1 * pkin(7) + t23; -t19, -t21, 0, t12; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
