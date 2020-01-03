% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4PRPR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:32
% EndTime: 2019-12-31 16:25:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (48->33), mult. (63->32), div. (0->0), fcn. (96->6), ass. (0->21)
t12 = sin(pkin(6));
t15 = sin(qJ(2));
t23 = qJ(3) * t15;
t17 = cos(qJ(2));
t5 = t12 * t17;
t28 = pkin(2) * t5 + t12 * t23;
t27 = t12 * t15;
t13 = cos(pkin(6));
t26 = t13 * t15;
t6 = t13 * t17;
t14 = sin(qJ(4));
t25 = t14 * t15;
t16 = cos(qJ(4));
t24 = t15 * t16;
t22 = t12 * pkin(1) + 0;
t11 = qJ(1) + 0;
t21 = t13 * pkin(1) + t12 * pkin(4) + 0;
t20 = pkin(2) * t6 + t13 * t23 + t21;
t19 = -t13 * pkin(4) + t22;
t18 = t15 * pkin(2) - t17 * qJ(3) + t11;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t13, -t12, 0, 0; t12, t13, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t6, -t26, t12, t21; t5, -t27, -t13, t19; t15, t17, 0, t11; 0, 0, 0, 1; t12, -t6, t26, t20; -t13, -t5, t27, t19 + t28; 0, -t15, -t17, t18; 0, 0, 0, 1; t12 * t16 + t13 * t25, -t12 * t14 + t13 * t24, t6, t12 * pkin(3) + pkin(5) * t6 + t20; t12 * t25 - t13 * t16, t12 * t24 + t13 * t14, t5, pkin(5) * t5 + (-pkin(3) - pkin(4)) * t13 + t22 + t28; -t17 * t14, -t17 * t16, t15, t15 * pkin(5) + t18; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
