% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRRP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:00
% EndTime: 2019-12-31 21:56:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (118->41), mult. (105->41), div. (0->0), fcn. (155->8), ass. (0->33)
t22 = qJ(2) + qJ(3);
t18 = sin(t22);
t23 = sin(qJ(4));
t42 = t18 * t23;
t25 = sin(qJ(1));
t12 = t25 * t18;
t19 = cos(t22);
t41 = t25 * t19;
t40 = t25 * t23;
t26 = cos(qJ(4));
t39 = t25 * t26;
t28 = cos(qJ(1));
t13 = t28 * t18;
t38 = t28 * t19;
t37 = t28 * t23;
t36 = t28 * t26;
t21 = pkin(5) + 0;
t24 = sin(qJ(2));
t35 = t24 * pkin(2) + t21;
t27 = cos(qJ(2));
t16 = t27 * pkin(2) + pkin(1);
t29 = -pkin(7) - pkin(6);
t34 = t25 * t16 + t28 * t29 + 0;
t33 = pkin(3) * t41 + pkin(8) * t12 + t34;
t32 = t28 * t16 - t25 * t29 + 0;
t31 = t18 * pkin(3) - t19 * pkin(8) + t35;
t30 = pkin(3) * t38 + pkin(8) * t13 + t32;
t11 = t18 * t26;
t4 = t19 * t36 + t40;
t3 = t19 * t37 - t39;
t2 = t19 * t39 - t37;
t1 = t19 * t40 + t36;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t28 * t27, -t28 * t24, t25, t28 * pkin(1) + t25 * pkin(6) + 0; t25 * t27, -t25 * t24, -t28, t25 * pkin(1) - t28 * pkin(6) + 0; t24, t27, 0, t21; 0, 0, 0, 1; t38, -t13, t25, t32; t41, -t12, -t28, t34; t18, t19, 0, t35; 0, 0, 0, 1; t4, -t3, t13, t30; t2, -t1, t12, t33; t11, -t42, -t19, t31; 0, 0, 0, 1; t4, t13, t3, t4 * pkin(4) + t3 * qJ(5) + t30; t2, t12, t1, t2 * pkin(4) + t1 * qJ(5) + t33; t11, -t19, t42, (pkin(4) * t26 + qJ(5) * t23) * t18 + t31; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
