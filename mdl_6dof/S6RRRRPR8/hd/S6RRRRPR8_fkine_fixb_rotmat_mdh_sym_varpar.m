% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:17:03
% EndTime: 2018-11-23 18:17:04
% DurationCPUTime: 0.16s
% Computational Cost: add. (204->64), mult. (230->71), div. (0->0), fcn. (314->10), ass. (0->39)
t25 = qJ(3) + qJ(4);
t19 = sin(t25);
t20 = cos(t25);
t33 = cos(qJ(1));
t29 = sin(qJ(1));
t32 = cos(qJ(2));
t49 = t29 * t32;
t3 = t19 * t49 + t33 * t20;
t4 = -t33 * t19 + t20 * t49;
t53 = t4 * pkin(4) + t3 * qJ(5);
t47 = t33 * t32;
t5 = t19 * t47 - t29 * t20;
t6 = t29 * t19 + t20 * t47;
t52 = t6 * pkin(4) + t5 * qJ(5);
t28 = sin(qJ(2));
t51 = t28 * t19;
t12 = t28 * t20;
t27 = sin(qJ(3));
t50 = t29 * t27;
t15 = t29 * t28;
t16 = t33 * t28;
t24 = pkin(6) + 0;
t45 = t29 * pkin(1) + 0;
t34 = -pkin(9) - pkin(8);
t44 = t28 * (-pkin(10) - t34);
t43 = t33 * pkin(1) + t29 * pkin(7) + 0;
t31 = cos(qJ(3));
t17 = t31 * pkin(3) + pkin(2);
t42 = t28 * t17 + t32 * t34 + t24;
t41 = pkin(2) * t32 + pkin(8) * t28;
t40 = -t33 * pkin(7) + t45;
t39 = pkin(3) * t50 + t17 * t47 + t43;
t38 = pkin(4) * t12 + qJ(5) * t51 + t42;
t37 = t17 * t49 + (-pkin(3) * t27 - pkin(7)) * t33 + t45;
t36 = -t34 * t16 + t39;
t35 = -t34 * t15 + t37;
t30 = cos(qJ(6));
t26 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t29, 0, 0; t29, t33, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t47, -t16, t29, t43; t49, -t15, -t33, t40; t28, t32, 0, t24; 0, 0, 0, 1; t31 * t47 + t50, -t27 * t47 + t29 * t31, t16, t41 * t33 + t43; -t33 * t27 + t31 * t49, -t27 * t49 - t33 * t31, t15, t41 * t29 + t40; t28 * t31, -t28 * t27, -t32, t28 * pkin(2) - t32 * pkin(8) + t24; 0, 0, 0, 1; t6, -t5, t16, t36; t4, -t3, t15, t35; t12, -t51, -t32, t42; 0, 0, 0, 1; t6, t16, t5, t36 + t52; t4, t15, t3, t35 + t53; t12, -t32, t51, t38; 0, 0, 0, 1; t5 * t26 + t6 * t30, -t6 * t26 + t5 * t30, -t16, t6 * pkin(5) + t33 * t44 + t39 + t52; t3 * t26 + t4 * t30, -t4 * t26 + t3 * t30, -t15, t4 * pkin(5) + t29 * t44 + t37 + t53; (t19 * t26 + t20 * t30) * t28 (t19 * t30 - t20 * t26) * t28, t32, pkin(5) * t12 + t32 * pkin(10) + t38; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
