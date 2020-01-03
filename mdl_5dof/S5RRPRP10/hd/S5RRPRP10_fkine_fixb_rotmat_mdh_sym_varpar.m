% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRP10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:16
% EndTime: 2019-12-31 20:09:17
% DurationCPUTime: 0.11s
% Computational Cost: add. (78->46), mult. (110->40), div. (0->0), fcn. (156->6), ass. (0->33)
t20 = sin(qJ(1));
t19 = sin(qJ(2));
t32 = qJ(3) * t19;
t22 = cos(qJ(2));
t9 = t20 * t22;
t40 = pkin(2) * t9 + t20 * t32;
t18 = sin(qJ(4));
t39 = pkin(4) * t18;
t38 = t20 * t19;
t21 = cos(qJ(4));
t37 = t20 * t21;
t36 = t22 * t18;
t35 = t22 * t21;
t23 = cos(qJ(1));
t34 = t23 * t19;
t33 = t23 * t21;
t10 = t23 * t22;
t16 = pkin(5) + 0;
t31 = t20 * pkin(1) + 0;
t30 = t19 * pkin(2) + t16;
t29 = t23 * pkin(1) + t20 * pkin(6) + 0;
t28 = t31 + t40;
t27 = -t23 * pkin(6) + t31;
t26 = pkin(2) * t10 + t23 * t32 + t29;
t17 = -qJ(5) - pkin(7);
t25 = -t17 * t22 + t19 * t39;
t24 = -t22 * qJ(3) + t30;
t11 = t21 * pkin(4) + pkin(3);
t4 = t18 * t38 - t33;
t3 = t23 * t18 + t19 * t37;
t2 = t18 * t34 + t37;
t1 = -t20 * t18 + t19 * t33;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t20, 0, 0; t20, t23, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t10, -t34, t20, t29; t9, -t38, -t23, t27; t19, t22, 0, t16; 0, 0, 0, 1; t20, -t10, t34, t26; -t23, -t9, t38, t27 + t40; 0, -t19, -t22, t24; 0, 0, 0, 1; t2, t1, t10, t20 * pkin(3) + pkin(7) * t10 + t26; t4, t3, t9, pkin(7) * t9 + (-pkin(3) - pkin(6)) * t23 + t28; -t36, -t35, t19, t19 * pkin(7) + t24; 0, 0, 0, 1; t2, t1, t10, t20 * t11 + t25 * t23 + t26; t4, t3, t9, (-pkin(6) - t11) * t23 + t25 * t20 + t28; -t36, -t35, t19, -t19 * t17 + (-qJ(3) - t39) * t22 + t30; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
