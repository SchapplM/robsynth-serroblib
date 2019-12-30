% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-29 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRPR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 16:56:04
% EndTime: 2019-12-29 16:56:04
% DurationCPUTime: 0.19s
% Computational Cost: add. (96->35), mult. (78->28), div. (0->0), fcn. (120->8), ass. (0->24)
t18 = sin(qJ(3));
t19 = sin(qJ(1));
t32 = t19 * t18;
t16 = pkin(5) + 0;
t31 = qJ(3) + pkin(8);
t30 = t19 * pkin(1) + 0;
t29 = -pkin(6) + t16;
t22 = cos(qJ(1));
t28 = t22 * pkin(1) + t19 * qJ(2) + 0;
t27 = cos(t31);
t26 = sin(t31);
t21 = cos(qJ(3));
t11 = t21 * pkin(3) + pkin(2);
t25 = pkin(3) * t32 + t22 * t11 + t28;
t24 = -t22 * qJ(2) + t30;
t23 = t19 * t11 + (-pkin(3) * t18 - qJ(2)) * t22 + t30;
t20 = cos(qJ(5));
t17 = sin(qJ(5));
t12 = -qJ(4) + t29;
t4 = -t22 * t18 + t19 * t21;
t3 = -t22 * t21 - t32;
t2 = -t19 * t27 + t22 * t26;
t1 = -t19 * t26 - t22 * t27;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t19, 0, 0; t19, t22, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t22, 0, t19, t28; t19, 0, -t22, t24; 0, 1, 0, t16; 0, 0, 0, 1; -t3, t4, 0, t22 * pkin(2) + t28; t4, t3, 0, t19 * pkin(2) + t24; 0, 0, -1, t29; 0, 0, 0, 1; -t1, -t2, 0, t25; -t2, t1, 0, t23; 0, 0, -1, t12; 0, 0, 0, 1; -t1 * t20, t1 * t17, t2, -t1 * pkin(4) + t2 * pkin(7) + t25; -t2 * t20, t2 * t17, -t1, -t2 * pkin(4) - t1 * pkin(7) + t23; -t17, -t20, 0, t12; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
