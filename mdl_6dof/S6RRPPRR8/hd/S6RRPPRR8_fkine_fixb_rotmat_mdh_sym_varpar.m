% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2018-11-23 16:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPPRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:53:53
% EndTime: 2018-11-23 16:53:53
% DurationCPUTime: 0.15s
% Computational Cost: add. (166->68), mult. (286->77), div. (0->0), fcn. (388->10), ass. (0->40)
t27 = cos(pkin(10));
t29 = sin(qJ(2));
t13 = t29 * t27;
t26 = sin(pkin(10));
t50 = t29 * t26;
t52 = pkin(3) * t13 + qJ(4) * t50;
t28 = sin(qJ(5));
t51 = t26 * t28;
t30 = sin(qJ(1));
t15 = t30 * t29;
t32 = cos(qJ(2));
t49 = t30 * t32;
t33 = cos(qJ(1));
t16 = t33 * t29;
t48 = t33 * t32;
t47 = qJ(3) * t29;
t24 = pkin(6) + 0;
t46 = t29 * pkin(2) + t24;
t45 = pkin(5) * t28 + qJ(4);
t44 = t33 * pkin(1) + t30 * pkin(7) + 0;
t43 = t46 + t52;
t42 = t30 * pkin(1) - t33 * pkin(7) + 0;
t41 = pkin(2) * t48 + t33 * t47 + t44;
t40 = -t32 * qJ(3) + t46;
t6 = t30 * t26 + t27 * t48;
t39 = t6 * pkin(3) + t41;
t38 = pkin(2) * t49 + t30 * t47 + t42;
t4 = -t33 * t26 + t27 * t49;
t37 = t4 * pkin(3) + t38;
t5 = t26 * t48 - t30 * t27;
t36 = t5 * qJ(4) + t39;
t3 = t26 * t49 + t33 * t27;
t35 = t3 * qJ(4) + t37;
t34 = -pkin(9) - pkin(8);
t31 = cos(qJ(5));
t25 = qJ(5) + qJ(6);
t19 = cos(t25);
t18 = sin(t25);
t17 = t31 * pkin(5) + pkin(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t30, 0, 0; t30, t33, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t48, -t16, t30, t44; t49, -t15, -t33, t42; t29, t32, 0, t24; 0, 0, 0, 1; t6, -t5, t16, t41; t4, -t3, t15, t38; t13, -t50, -t32, t40; 0, 0, 0, 1; t6, t16, t5, t36; t4, t15, t3, t35; t13, -t32, t50, t40 + t52; 0, 0, 0, 1; t5 * t28 + t6 * t31, -t6 * t28 + t5 * t31, -t16, t6 * pkin(4) - pkin(8) * t16 + t36; t3 * t28 + t4 * t31, -t4 * t28 + t3 * t31, -t15, t4 * pkin(4) - pkin(8) * t15 + t35; (t27 * t31 + t51) * t29 (t26 * t31 - t27 * t28) * t29, t32, pkin(4) * t13 + (pkin(8) - qJ(3)) * t32 + t43; 0, 0, 0, 1; t5 * t18 + t6 * t19, -t6 * t18 + t5 * t19, -t16, t34 * t16 + t6 * t17 + t45 * t5 + t39; t3 * t18 + t4 * t19, -t4 * t18 + t3 * t19, -t15, t34 * t15 + t4 * t17 + t45 * t3 + t37; (t18 * t26 + t19 * t27) * t29 (-t18 * t27 + t19 * t26) * t29, t32 (-qJ(3) - t34) * t32 + (pkin(5) * t51 + t17 * t27) * t29 + t43; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
