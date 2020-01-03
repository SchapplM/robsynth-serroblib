% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4PRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:08
% EndTime: 2019-12-31 16:30:08
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->30), mult. (88->33), div. (0->0), fcn. (130->6), ass. (0->25)
t18 = sin(pkin(6));
t21 = sin(qJ(2));
t9 = t18 * t21;
t23 = cos(qJ(2));
t33 = t18 * t23;
t19 = cos(pkin(6));
t11 = t19 * t21;
t32 = t19 * t23;
t20 = sin(qJ(3));
t31 = t20 * t23;
t30 = t21 * t20;
t22 = cos(qJ(3));
t29 = t22 * t23;
t17 = qJ(1) + 0;
t28 = t19 * pkin(1) + t18 * pkin(4) + 0;
t27 = t18 * pkin(1) - t19 * pkin(4) + 0;
t26 = pkin(2) * t32 + pkin(5) * t11 + t28;
t25 = t21 * pkin(2) - t23 * pkin(5) + t17;
t24 = pkin(2) * t33 + pkin(5) * t9 + t27;
t12 = t21 * t22;
t4 = t18 * t20 + t19 * t29;
t3 = -t18 * t22 + t19 * t31;
t2 = t18 * t29 - t19 * t20;
t1 = t18 * t31 + t19 * t22;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t18, 0, 0; t18, t19, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; t32, -t11, t18, t28; t33, -t9, -t19, t27; t21, t23, 0, t17; 0, 0, 0, 1; t4, -t3, t11, t26; t2, -t1, t9, t24; t12, -t30, -t23, t25; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(3) + t3 * qJ(4) + t26; t2, t9, t1, t2 * pkin(3) + t1 * qJ(4) + t24; t12, -t23, t30, (pkin(3) * t22 + qJ(4) * t20) * t21 + t25; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
