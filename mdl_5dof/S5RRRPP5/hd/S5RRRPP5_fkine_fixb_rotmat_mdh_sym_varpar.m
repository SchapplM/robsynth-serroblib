% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-29 19:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRPP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 19:44:19
% EndTime: 2019-12-29 19:44:19
% DurationCPUTime: 0.19s
% Computational Cost: add. (97->37), mult. (68->30), div. (0->0), fcn. (104->6), ass. (0->23)
t18 = qJ(2) + qJ(3);
t15 = cos(t18);
t22 = cos(qJ(1));
t10 = t22 * t15;
t14 = sin(t18);
t30 = qJ(4) * t14;
t31 = pkin(3) * t10 + t22 * t30;
t20 = sin(qJ(1));
t8 = t20 * t15;
t17 = pkin(5) + 0;
t21 = cos(qJ(2));
t12 = t21 * pkin(2) + pkin(1);
t29 = t22 * t12 + 0;
t19 = sin(qJ(2));
t28 = t19 * pkin(2) + t17;
t23 = -pkin(7) - pkin(6);
t27 = t20 * t12 + t22 * t23 + 0;
t26 = pkin(3) * t8 + t20 * t30 + t27;
t25 = -t20 * t23 + t29;
t24 = t14 * pkin(3) - t15 * qJ(4) + t28;
t9 = t22 * t14;
t7 = t20 * t14;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t20, 0, 0; t20, t22, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; t22 * t21, -t22 * t19, t20, t22 * pkin(1) + t20 * pkin(6) + 0; t20 * t21, -t20 * t19, -t22, t20 * pkin(1) - t22 * pkin(6) + 0; t19, t21, 0, t17; 0, 0, 0, 1; t10, -t9, t20, t25; t8, -t7, -t22, t27; t14, t15, 0, t28; 0, 0, 0, 1; t10, t20, t9, t25 + t31; t8, -t22, t7, t26; t14, 0, -t15, t24; 0, 0, 0, 1; t10, t9, -t20, pkin(4) * t10 + (-qJ(5) - t23) * t20 + t29 + t31; t8, t7, t22, pkin(4) * t8 + t22 * qJ(5) + t26; t14, -t15, 0, t14 * pkin(4) + t24; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
