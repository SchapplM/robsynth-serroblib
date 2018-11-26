% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:27:39
% EndTime: 2018-11-23 18:27:39
% DurationCPUTime: 0.14s
% Computational Cost: add. (189->61), mult. (137->64), div. (0->0), fcn. (196->10), ass. (0->43)
t31 = -pkin(10) - pkin(9);
t25 = sin(qJ(4));
t47 = t25 * pkin(4);
t28 = cos(qJ(4));
t12 = t28 * pkin(4) + pkin(3);
t23 = qJ(4) + qJ(5);
t15 = sin(t23);
t24 = qJ(2) + qJ(3);
t16 = sin(t24);
t46 = t16 * t15;
t18 = cos(t24);
t27 = sin(qJ(1));
t45 = t27 * t18;
t44 = t27 * t25;
t43 = t27 * t28;
t30 = cos(qJ(1));
t42 = t30 * t18;
t41 = t30 * t25;
t40 = t30 * t28;
t22 = pkin(6) + 0;
t29 = cos(qJ(2));
t13 = t29 * pkin(2) + pkin(1);
t39 = t30 * t13 + 0;
t26 = sin(qJ(2));
t38 = t26 * pkin(2) + t22;
t32 = -pkin(8) - pkin(7);
t37 = t27 * t13 + t30 * t32 + 0;
t36 = pkin(3) * t18 + pkin(9) * t16;
t21 = -qJ(6) + t31;
t17 = cos(t23);
t5 = pkin(5) * t17 + t12;
t35 = -t16 * t21 + t18 * t5;
t34 = t12 * t18 - t16 * t31;
t33 = -t27 * t32 + t39;
t11 = t30 * t16;
t10 = t27 * t16;
t7 = t16 * t17;
t6 = pkin(5) * t15 + t47;
t4 = t27 * t15 + t17 * t42;
t3 = -t15 * t42 + t27 * t17;
t2 = -t30 * t15 + t17 * t45;
t1 = -t15 * t45 - t30 * t17;
t8 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t27, 0, 0; t27, t30, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t30 * t29, -t30 * t26, t27, t30 * pkin(1) + t27 * pkin(7) + 0; t27 * t29, -t27 * t26, -t30, t27 * pkin(1) - t30 * pkin(7) + 0; t26, t29, 0, t22; 0, 0, 0, 1; t42, -t11, t27, t33; t45, -t10, -t30, t37; t16, t18, 0, t38; 0, 0, 0, 1; t18 * t40 + t44, -t18 * t41 + t43, t11, t36 * t30 + t33; t18 * t43 - t41, -t18 * t44 - t40, t10, t36 * t27 + t37; t16 * t28, -t16 * t25, -t18, t16 * pkin(3) - t18 * pkin(9) + t38; 0, 0, 0, 1; t4, t3, t11, t34 * t30 + (-t32 + t47) * t27 + t39; t2, t1, t10, -pkin(4) * t41 + t34 * t27 + t37; t7, -t46, -t18, t16 * t12 + t18 * t31 + t38; 0, 0, 0, 1; t4, t3, t11, t35 * t30 + (-t32 + t6) * t27 + t39; t2, t1, t10, t35 * t27 - t30 * t6 + t37; t7, -t46, -t18, t16 * t5 + t18 * t21 + t38; 0, 0, 0, 1;];
T_ges = t8;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
