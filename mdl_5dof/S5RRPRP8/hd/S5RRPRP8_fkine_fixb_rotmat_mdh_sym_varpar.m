% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP8
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRP8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:13
% EndTime: 2019-12-31 20:03:13
% DurationCPUTime: 0.13s
% Computational Cost: add. (77->41), mult. (119->38), div. (0->0), fcn. (167->6), ass. (0->31)
t22 = sin(qJ(1));
t24 = cos(qJ(2));
t11 = t22 * t24;
t21 = sin(qJ(2));
t35 = qJ(3) * t21;
t39 = pkin(2) * t11 + t22 * t35;
t20 = sin(qJ(4));
t38 = t21 * t20;
t37 = t22 * t21;
t25 = cos(qJ(1));
t36 = t25 * t21;
t12 = t25 * t24;
t18 = pkin(5) + 0;
t34 = t22 * pkin(1) + 0;
t33 = t21 * pkin(2) + t18;
t32 = t25 * pkin(1) + t22 * pkin(6) + 0;
t31 = t34 + t39;
t23 = cos(qJ(4));
t6 = -t24 * t20 + t21 * t23;
t30 = t24 * t23 + t38;
t29 = -t25 * pkin(6) + t34;
t28 = pkin(2) * t12 + t25 * t35 + t32;
t13 = t23 * pkin(4) + pkin(3);
t27 = pkin(4) * t38 + t13 * t24;
t26 = -t24 * qJ(3) + t33;
t19 = -qJ(5) - pkin(7);
t4 = t30 * t25;
t3 = t6 * t25;
t2 = t30 * t22;
t1 = t6 * t22;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, -t22, 0, 0; t22, t25, 0, 0; 0, 0, 1, t18; 0, 0, 0, 1; t12, -t36, t22, t32; t11, -t37, -t25, t29; t21, t24, 0, t18; 0, 0, 0, 1; t12, t22, t36, t28; t11, -t25, t37, t29 + t39; t21, 0, -t24, t26; 0, 0, 0, 1; t4, t3, -t22, pkin(3) * t12 - t22 * pkin(7) + t28; t2, t1, t25, pkin(3) * t11 + (-pkin(6) + pkin(7)) * t25 + t31; t6, -t30, 0, t21 * pkin(3) + t26; 0, 0, 0, 1; t4, t3, -t22, t22 * t19 + t25 * t27 + t28; t2, t1, t25, (-pkin(6) - t19) * t25 + t27 * t22 + t31; t6, -t30, 0, t21 * t13 + (-pkin(4) * t20 - qJ(3)) * t24 + t33; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
