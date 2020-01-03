% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPPR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:17
% EndTime: 2019-12-31 19:46:17
% DurationCPUTime: 0.11s
% Computational Cost: add. (88->51), mult. (110->51), div. (0->0), fcn. (156->8), ass. (0->29)
t20 = sin(qJ(1));
t19 = sin(qJ(2));
t32 = qJ(3) * t19;
t21 = cos(qJ(2));
t6 = t20 * t21;
t36 = pkin(2) * t6 + t20 * t32;
t16 = sin(pkin(8));
t35 = pkin(4) * t16;
t34 = t20 * t19;
t22 = cos(qJ(1));
t33 = t22 * t19;
t7 = t22 * t21;
t31 = qJ(4) * t21;
t15 = pkin(5) + 0;
t30 = t20 * pkin(1) + 0;
t29 = t19 * pkin(2) + t15;
t28 = t22 * pkin(1) + t20 * pkin(6) + 0;
t27 = t30 + t36;
t26 = -t22 * pkin(6) + t30;
t25 = pkin(2) * t7 + t22 * t32 + t28;
t18 = -pkin(7) - qJ(4);
t24 = -t18 * t21 + t19 * t35;
t23 = -t21 * qJ(3) + t29;
t17 = cos(pkin(8));
t14 = pkin(8) + qJ(5);
t9 = cos(t14);
t8 = sin(t14);
t5 = t17 * pkin(4) + pkin(3);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t20, 0, 0; t20, t22, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t7, -t33, t20, t28; t6, -t34, -t22, t26; t19, t21, 0, t15; 0, 0, 0, 1; t20, -t7, t33, t25; -t22, -t6, t34, t26 + t36; 0, -t19, -t21, t23; 0, 0, 0, 1; t16 * t33 + t20 * t17, -t20 * t16 + t17 * t33, t7, t20 * pkin(3) + t22 * t31 + t25; t16 * t34 - t22 * t17, t22 * t16 + t17 * t34, t6, t20 * t31 + (-pkin(3) - pkin(6)) * t22 + t27; -t21 * t16, -t21 * t17, t19, t19 * qJ(4) + t23; 0, 0, 0, 1; t20 * t9 + t8 * t33, -t20 * t8 + t9 * t33, t7, t20 * t5 + t24 * t22 + t25; -t22 * t9 + t8 * t34, t22 * t8 + t9 * t34, t6, (-pkin(6) - t5) * t22 + t24 * t20 + t27; -t21 * t8, -t21 * t9, t19, -t19 * t18 + (-qJ(3) - t35) * t21 + t29; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
