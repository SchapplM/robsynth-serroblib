% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-29 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRP11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 17:31:21
% EndTime: 2019-12-29 17:31:22
% DurationCPUTime: 0.26s
% Computational Cost: add. (118->41), mult. (105->41), div. (0->0), fcn. (155->8), ass. (0->33)
t21 = pkin(8) + qJ(3);
t18 = sin(t21);
t26 = sin(qJ(4));
t42 = t18 * t26;
t27 = sin(qJ(1));
t12 = t27 * t18;
t19 = cos(t21);
t41 = t27 * t19;
t40 = t27 * t26;
t28 = cos(qJ(4));
t39 = t27 * t28;
t29 = cos(qJ(1));
t13 = t29 * t18;
t38 = t29 * t19;
t37 = t29 * t26;
t36 = t29 * t28;
t22 = pkin(5) + 0;
t23 = sin(pkin(8));
t35 = t23 * pkin(2) + t22;
t24 = cos(pkin(8));
t15 = t24 * pkin(2) + pkin(1);
t25 = -pkin(6) - qJ(2);
t34 = t27 * t15 + t29 * t25 + 0;
t33 = pkin(3) * t41 + pkin(7) * t12 + t34;
t32 = t29 * t15 - t27 * t25 + 0;
t31 = t18 * pkin(3) - t19 * pkin(7) + t35;
t30 = pkin(3) * t38 + pkin(7) * t13 + t32;
t11 = t18 * t28;
t4 = t19 * t36 + t40;
t3 = t19 * t37 - t39;
t2 = t19 * t39 - t37;
t1 = t19 * t40 + t36;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t27, 0, 0; t27, t29, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t29 * t24, -t29 * t23, t27, t29 * pkin(1) + t27 * qJ(2) + 0; t27 * t24, -t27 * t23, -t29, t27 * pkin(1) - t29 * qJ(2) + 0; t23, t24, 0, t22; 0, 0, 0, 1; t38, -t13, t27, t32; t41, -t12, -t29, t34; t18, t19, 0, t35; 0, 0, 0, 1; t4, -t3, t13, t30; t2, -t1, t12, t33; t11, -t42, -t19, t31; 0, 0, 0, 1; t4, t13, t3, t4 * pkin(4) + t3 * qJ(5) + t30; t2, t12, t1, t2 * pkin(4) + t1 * qJ(5) + t33; t11, -t19, t42, (pkin(4) * t28 + qJ(5) * t26) * t18 + t31; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
