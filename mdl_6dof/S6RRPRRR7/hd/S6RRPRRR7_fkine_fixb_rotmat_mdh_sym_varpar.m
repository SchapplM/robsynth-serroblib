% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:25:01
% EndTime: 2018-11-23 17:25:01
% DurationCPUTime: 0.15s
% Computational Cost: add. (148->61), mult. (248->66), div. (0->0), fcn. (338->10), ass. (0->40)
t31 = sin(qJ(1));
t34 = cos(qJ(2));
t15 = t31 * t34;
t30 = sin(qJ(2));
t47 = qJ(3) * t30;
t53 = pkin(2) * t15 + t31 * t47;
t35 = cos(qJ(1));
t52 = pkin(3) * t15 + pkin(8) * t35;
t28 = sin(qJ(5));
t51 = pkin(5) * t28;
t50 = t31 * t30;
t29 = sin(qJ(4));
t49 = t34 * t29;
t48 = t35 * t30;
t16 = t35 * t34;
t26 = pkin(6) + 0;
t46 = pkin(1) * t31 + 0;
t45 = pkin(1) * t35 + pkin(7) * t31 + 0;
t33 = cos(qJ(4));
t5 = t29 * t30 + t33 * t34;
t44 = -pkin(7) * t35 + t46;
t43 = pkin(2) * t16 + t35 * t47 + t45;
t42 = pkin(2) * t30 - qJ(3) * t34 + t26;
t41 = pkin(3) * t16 + t43;
t40 = t44 + t53;
t39 = pkin(3) * t30 + t42;
t38 = t40 + t52;
t37 = -pkin(8) * t31 + t41;
t36 = -pkin(10) - pkin(9);
t32 = cos(qJ(5));
t27 = qJ(5) + qJ(6);
t19 = cos(t27);
t18 = sin(t27);
t17 = pkin(5) * t32 + pkin(4);
t6 = t30 * t33 - t49;
t4 = t5 * t35;
t3 = t16 * t29 - t33 * t48;
t2 = t5 * t31;
t1 = t31 * t49 - t33 * t50;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t31, 0, 0; t31, t35, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t16, -t48, t31, t45; t15, -t50, -t35, t44; t30, t34, 0, t26; 0, 0, 0, 1; t16, t31, t48, t43; t15, -t35, t50, t40; t30, 0, -t34, t42; 0, 0, 0, 1; t4, -t3, -t31, t37; t2, -t1, t35, t38; t6, -t5, 0, t39; 0, 0, 0, 1; -t28 * t31 + t32 * t4, -t28 * t4 - t31 * t32, t3, pkin(4) * t4 + pkin(9) * t3 + t37; t2 * t32 + t28 * t35, -t2 * t28 + t32 * t35, t1, pkin(4) * t2 + pkin(9) * t1 + t38; t6 * t32, -t6 * t28, t5, pkin(4) * t6 + pkin(9) * t5 + t39; 0, 0, 0, 1; -t18 * t31 + t19 * t4, -t18 * t4 - t19 * t31, t3, t4 * t17 - t3 * t36 + (-pkin(8) - t51) * t31 + t41; t18 * t35 + t19 * t2, -t18 * t2 + t19 * t35, t1, -t1 * t36 + t2 * t17 + (-pkin(7) + t51) * t35 + t46 + t52 + t53; t6 * t19, -t6 * t18, t5, t17 * t6 - t36 * t5 + t39; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
