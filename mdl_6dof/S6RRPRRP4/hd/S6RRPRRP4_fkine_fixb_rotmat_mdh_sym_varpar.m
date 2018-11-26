% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:12:17
% EndTime: 2018-11-23 17:12:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (201->57), mult. (152->60), div. (0->0), fcn. (215->10), ass. (0->43)
t24 = qJ(2) + pkin(10);
t19 = sin(t24);
t26 = qJ(4) + qJ(5);
t21 = sin(t26);
t52 = t19 * t21;
t30 = sin(qJ(1));
t11 = t30 * t19;
t20 = cos(t24);
t51 = t30 * t20;
t50 = t30 * t21;
t22 = cos(t26);
t49 = t30 * t22;
t28 = sin(qJ(4));
t48 = t30 * t28;
t31 = cos(qJ(4));
t47 = t30 * t31;
t33 = cos(qJ(1));
t12 = t33 * t19;
t46 = t33 * t20;
t45 = t33 * t21;
t44 = t33 * t22;
t43 = t33 * t28;
t42 = t33 * t31;
t25 = pkin(6) + 0;
t29 = sin(qJ(2));
t41 = t29 * pkin(2) + t25;
t32 = cos(qJ(2));
t18 = t32 * pkin(2) + pkin(1);
t27 = -qJ(3) - pkin(7);
t40 = t30 * t18 + t33 * t27 + 0;
t39 = pkin(3) * t20 + pkin(8) * t19;
t17 = t31 * pkin(4) + pkin(3);
t34 = -pkin(9) - pkin(8);
t38 = t19 * t17 + t20 * t34 + t41;
t37 = t33 * t18 - t30 * t27 + 0;
t36 = pkin(4) * t48 - t34 * t12 + t17 * t46 + t37;
t35 = -pkin(4) * t43 - t34 * t11 + t17 * t51 + t40;
t8 = t19 * t22;
t4 = t20 * t44 + t50;
t3 = t20 * t45 - t49;
t2 = t20 * t49 - t45;
t1 = t20 * t50 + t44;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t30, 0, 0; t30, t33, 0, 0; 0, 0, 1, t25; 0, 0, 0, 1; t33 * t32, -t33 * t29, t30, t33 * pkin(1) + t30 * pkin(7) + 0; t30 * t32, -t30 * t29, -t33, t30 * pkin(1) - t33 * pkin(7) + 0; t29, t32, 0, t25; 0, 0, 0, 1; t46, -t12, t30, t37; t51, -t11, -t33, t40; t19, t20, 0, t41; 0, 0, 0, 1; t20 * t42 + t48, -t20 * t43 + t47, t12, t39 * t33 + t37; t20 * t47 - t43, -t20 * t48 - t42, t11, t39 * t30 + t40; t19 * t31, -t19 * t28, -t20, t19 * pkin(3) - t20 * pkin(8) + t41; 0, 0, 0, 1; t4, -t3, t12, t36; t2, -t1, t11, t35; t8, -t52, -t20, t38; 0, 0, 0, 1; t4, t12, t3, t4 * pkin(5) + t3 * qJ(6) + t36; t2, t11, t1, t2 * pkin(5) + t1 * qJ(6) + t35; t8, -t20, t52 (pkin(5) * t22 + qJ(6) * t21) * t19 + t38; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
