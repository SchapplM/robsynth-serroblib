% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:25:05
% EndTime: 2018-11-23 16:25:05
% DurationCPUTime: 0.12s
% Computational Cost: add. (213->52), mult. (142->58), div. (0->0), fcn. (201->10), ass. (0->42)
t25 = qJ(1) + pkin(10);
t18 = sin(t25);
t27 = sin(qJ(4));
t52 = t18 * t27;
t31 = cos(qJ(3));
t51 = t18 * t31;
t19 = cos(t25);
t50 = t19 * t31;
t26 = qJ(4) + qJ(5);
t21 = sin(t26);
t49 = t21 * t31;
t22 = cos(t26);
t48 = t22 * t31;
t47 = t27 * t31;
t28 = sin(qJ(3));
t46 = t28 * t21;
t33 = -pkin(9) - pkin(8);
t45 = t28 * t33;
t30 = cos(qJ(4));
t44 = t30 * t31;
t43 = pkin(6) + 0;
t29 = sin(qJ(1));
t42 = t29 * pkin(1) + 0;
t32 = cos(qJ(1));
t41 = t32 * pkin(1) + 0;
t20 = qJ(2) + t43;
t40 = t18 * pkin(2) + t42;
t39 = t19 * pkin(2) + t18 * pkin(7) + t41;
t38 = pkin(3) * t31 + pkin(8) * t28;
t16 = t30 * pkin(4) + pkin(3);
t37 = t28 * t16 + t31 * t33 + t20;
t36 = -t19 * pkin(7) + t40;
t35 = pkin(4) * t52 + t16 * t50 - t19 * t45 + t39;
t34 = -t18 * t45 + t16 * t51 + (-pkin(4) * t27 - pkin(7)) * t19 + t40;
t12 = t28 * t22;
t11 = t19 * t28;
t10 = t18 * t28;
t4 = t18 * t21 + t19 * t48;
t3 = -t18 * t22 + t19 * t49;
t2 = t18 * t48 - t19 * t21;
t1 = t18 * t49 + t19 * t22;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t32, -t29, 0, 0; t29, t32, 0, 0; 0, 0, 1, t43; 0, 0, 0, 1; t19, -t18, 0, t41; t18, t19, 0, t42; 0, 0, 1, t20; 0, 0, 0, 1; t50, -t11, t18, t39; t51, -t10, -t19, t36; t28, t31, 0, t20; 0, 0, 0, 1; t19 * t44 + t52, t18 * t30 - t19 * t47, t11, t38 * t19 + t39; t18 * t44 - t19 * t27, -t18 * t47 - t19 * t30, t10, t38 * t18 + t36; t28 * t30, -t28 * t27, -t31, t28 * pkin(3) - t31 * pkin(8) + t20; 0, 0, 0, 1; t4, -t3, t11, t35; t2, -t1, t10, t34; t12, -t46, -t31, t37; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(5) + t3 * qJ(6) + t35; t2, t10, t1, t2 * pkin(5) + t1 * qJ(6) + t34; t12, -t31, t46 (pkin(5) * t22 + qJ(6) * t21) * t28 + t37; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
