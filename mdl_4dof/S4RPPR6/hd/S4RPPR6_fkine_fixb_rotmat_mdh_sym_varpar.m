% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-29 12:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RPPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 12:43:55
% EndTime: 2019-12-29 12:43:55
% DurationCPUTime: 0.16s
% Computational Cost: add. (46->28), mult. (67->28), div. (0->0), fcn. (101->6), ass. (0->21)
t15 = sin(qJ(1));
t12 = sin(pkin(6));
t25 = qJ(3) * t12;
t13 = cos(pkin(6));
t5 = t15 * t13;
t28 = pkin(2) * t5 + t15 * t25;
t27 = t15 * t12;
t17 = cos(qJ(1));
t26 = t17 * t12;
t6 = t17 * t13;
t11 = pkin(4) + 0;
t24 = t15 * pkin(1) + 0;
t23 = t17 * pkin(1) + t15 * qJ(2) + 0;
t22 = pkin(2) * t6 + t17 * t25 + t23;
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t21 = t12 * t16 - t13 * t14;
t20 = t12 * t14 + t13 * t16;
t19 = -t17 * qJ(2) + t24;
t18 = t12 * pkin(2) - t13 * qJ(3) + t11;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t15, 0, 0; t15, t17, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t6, -t26, t15, t23; t5, -t27, -t17, t19; t12, t13, 0, t11; 0, 0, 0, 1; t6, t15, t26, t22; t5, -t17, t27, t19 + t28; t12, 0, -t13, t18; 0, 0, 0, 1; t20 * t17, t21 * t17, -t15, pkin(3) * t6 - t15 * pkin(5) + t22; t20 * t15, t21 * t15, t17, pkin(3) * t5 + (pkin(5) - qJ(2)) * t17 + t24 + t28; t21, -t20, 0, t12 * pkin(3) + t18; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
