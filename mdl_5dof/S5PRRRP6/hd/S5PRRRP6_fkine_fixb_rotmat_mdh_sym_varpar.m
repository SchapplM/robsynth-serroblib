% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP6
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:38
% EndTime: 2019-12-05 16:50:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (116->47), mult. (132->54), div. (0->0), fcn. (187->8), ass. (0->34)
t22 = sin(pkin(8));
t24 = sin(qJ(3));
t42 = t22 * t24;
t27 = cos(qJ(2));
t41 = t22 * t27;
t23 = cos(pkin(8));
t40 = t23 * t27;
t39 = t24 * t27;
t21 = qJ(3) + qJ(4);
t15 = sin(t21);
t25 = sin(qJ(2));
t38 = t25 * t15;
t28 = -pkin(7) - pkin(6);
t37 = t25 * t28;
t26 = cos(qJ(3));
t36 = t26 * t27;
t35 = t22 * pkin(1) + 0;
t20 = qJ(1) + 0;
t34 = t23 * pkin(1) + t22 * pkin(5) + 0;
t13 = t26 * pkin(3) + pkin(2);
t33 = t25 * t13 + t27 * t28 + t20;
t32 = pkin(2) * t27 + pkin(6) * t25;
t31 = -t23 * pkin(5) + t35;
t30 = pkin(3) * t42 + t13 * t40 - t23 * t37 + t34;
t29 = -t22 * t37 + t13 * t41 + (-pkin(3) * t24 - pkin(5)) * t23 + t35;
t16 = cos(t21);
t12 = t23 * t25;
t11 = t22 * t25;
t9 = t25 * t16;
t4 = t22 * t15 + t16 * t40;
t3 = t15 * t40 - t22 * t16;
t2 = -t23 * t15 + t16 * t41;
t1 = t15 * t41 + t23 * t16;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t22, 0, 0; t22, t23, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t40, -t12, t22, t34; t41, -t11, -t23, t31; t25, t27, 0, t20; 0, 0, 0, 1; t23 * t36 + t42, t22 * t26 - t23 * t39, t12, t32 * t23 + t34; t22 * t36 - t23 * t24, -t22 * t39 - t23 * t26, t11, t32 * t22 + t31; t25 * t26, -t25 * t24, -t27, t25 * pkin(2) - t27 * pkin(6) + t20; 0, 0, 0, 1; t4, -t3, t12, t30; t2, -t1, t11, t29; t9, -t38, -t27, t33; 0, 0, 0, 1; t4, t12, t3, t4 * pkin(4) + t3 * qJ(5) + t30; t2, t11, t1, t2 * pkin(4) + t1 * qJ(5) + t29; t9, -t27, t38, (pkin(4) * t16 + qJ(5) * t15) * t25 + t33; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
