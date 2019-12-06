% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:29
% EndTime: 2019-12-05 16:45:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (118->41), mult. (105->41), div. (0->0), fcn. (155->8), ass. (0->33)
t22 = qJ(2) + qJ(3);
t18 = sin(t22);
t25 = sin(qJ(4));
t42 = t18 * t25;
t23 = sin(pkin(8));
t11 = t23 * t18;
t19 = cos(t22);
t41 = t23 * t19;
t40 = t23 * t25;
t27 = cos(qJ(4));
t39 = t23 * t27;
t24 = cos(pkin(8));
t12 = t24 * t18;
t38 = t24 * t19;
t37 = t24 * t25;
t36 = t24 * t27;
t21 = qJ(1) + 0;
t28 = cos(qJ(2));
t17 = t28 * pkin(2) + pkin(1);
t29 = -pkin(6) - pkin(5);
t35 = t23 * t17 + t24 * t29 + 0;
t26 = sin(qJ(2));
t34 = t26 * pkin(2) + t21;
t33 = pkin(3) * t41 + pkin(7) * t11 + t35;
t32 = t24 * t17 - t23 * t29 + 0;
t31 = pkin(3) * t38 + pkin(7) * t12 + t32;
t30 = t18 * pkin(3) - t19 * pkin(7) + t34;
t13 = t18 * t27;
t4 = t19 * t36 + t40;
t3 = t19 * t37 - t39;
t2 = t19 * t39 - t37;
t1 = t19 * t40 + t36;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t23, 0, 0; t23, t24, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t24 * t28, -t24 * t26, t23, t24 * pkin(1) + t23 * pkin(5) + 0; t23 * t28, -t23 * t26, -t24, t23 * pkin(1) - t24 * pkin(5) + 0; t26, t28, 0, t21; 0, 0, 0, 1; t38, -t12, t23, t32; t41, -t11, -t24, t35; t18, t19, 0, t34; 0, 0, 0, 1; t4, -t3, t12, t31; t2, -t1, t11, t33; t13, -t42, -t19, t30; 0, 0, 0, 1; t4, t12, t3, t4 * pkin(4) + t3 * qJ(5) + t31; t2, t11, t1, t2 * pkin(4) + t1 * qJ(5) + t33; t13, -t19, t42, (pkin(4) * t27 + qJ(5) * t25) * t18 + t30; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
