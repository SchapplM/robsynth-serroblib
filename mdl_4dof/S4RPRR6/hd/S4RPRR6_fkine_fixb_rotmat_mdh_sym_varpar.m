% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RPRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 13:11:48
% EndTime: 2019-12-29 13:11:48
% DurationCPUTime: 0.15s
% Computational Cost: add. (62->31), mult. (33->28), div. (0->0), fcn. (61->8), ass. (0->17)
t14 = cos(pkin(7));
t4 = t14 * pkin(2) + pkin(1);
t15 = -pkin(5) - qJ(2);
t12 = pkin(4) + 0;
t11 = pkin(7) + qJ(3);
t13 = sin(pkin(7));
t18 = t13 * pkin(2) + t12;
t17 = cos(qJ(1));
t16 = sin(qJ(1));
t10 = -pkin(6) + t15;
t7 = qJ(4) + t11;
t6 = cos(t11);
t5 = sin(t11);
t3 = cos(t7);
t2 = sin(t7);
t1 = pkin(3) * t6 + t4;
t8 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t16, 0, 0; t16, t17, 0, 0; 0, 0, 1, t12; 0, 0, 0, 1; t17 * t14, -t17 * t13, t16, t17 * pkin(1) + t16 * qJ(2) + 0; t16 * t14, -t16 * t13, -t17, t16 * pkin(1) - t17 * qJ(2) + 0; t13, t14, 0, t12; 0, 0, 0, 1; t17 * t6, -t17 * t5, t16, -t16 * t15 + t17 * t4 + 0; t16 * t6, -t16 * t5, -t17, t17 * t15 + t16 * t4 + 0; t5, t6, 0, t18; 0, 0, 0, 1; t17 * t3, -t17 * t2, t16, t17 * t1 - t16 * t10 + 0; t16 * t3, -t16 * t2, -t17, t16 * t1 + t17 * t10 + 0; t2, t3, 0, pkin(3) * t5 + t18; 0, 0, 0, 1;];
T_ges = t8;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
