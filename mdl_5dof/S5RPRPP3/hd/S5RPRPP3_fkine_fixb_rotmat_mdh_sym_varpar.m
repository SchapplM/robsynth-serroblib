% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRPP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:13
% EndTime: 2019-12-31 18:12:14
% DurationCPUTime: 0.09s
% Computational Cost: add. (100->40), mult. (68->31), div. (0->0), fcn. (104->6), ass. (0->24)
t17 = pkin(7) + qJ(3);
t15 = cos(t17);
t23 = cos(qJ(1));
t10 = t23 * t15;
t14 = sin(t17);
t31 = qJ(4) * t14;
t32 = pkin(3) * t10 + t23 * t31;
t22 = sin(qJ(1));
t8 = t22 * t15;
t30 = qJ(5) * t15;
t18 = pkin(5) + 0;
t20 = cos(pkin(7));
t12 = t20 * pkin(2) + pkin(1);
t29 = t23 * t12 + 0;
t19 = sin(pkin(7));
t28 = t19 * pkin(2) + t18;
t21 = -pkin(6) - qJ(2);
t27 = t22 * t12 + t23 * t21 + 0;
t26 = pkin(3) * t8 + t22 * t31 + t27;
t25 = -t22 * t21 + t29;
t24 = t14 * pkin(3) - t15 * qJ(4) + t28;
t9 = t23 * t14;
t7 = t22 * t14;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t22, 0, 0; t22, t23, 0, 0; 0, 0, 1, t18; 0, 0, 0, 1; t23 * t20, -t23 * t19, t22, t23 * pkin(1) + t22 * qJ(2) + 0; t22 * t20, -t22 * t19, -t23, t22 * pkin(1) - t23 * qJ(2) + 0; t19, t20, 0, t18; 0, 0, 0, 1; t10, -t9, t22, t25; t8, -t7, -t23, t27; t14, t15, 0, t28; 0, 0, 0, 1; t22, -t10, t9, t25 + t32; -t23, -t8, t7, t26; 0, -t14, -t15, t24; 0, 0, 0, 1; t22, t9, t10, t23 * t30 + (pkin(4) - t21) * t22 + t29 + t32; -t23, t7, t8, -t23 * pkin(4) + t22 * t30 + t26; 0, -t15, t14, t14 * qJ(5) + t24; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
