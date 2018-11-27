% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% T_c_mdh [4x4x(7+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   8:  mdh base (link 0) -> mdh frame (8-1), link (8-1)
%   ...
%   7+1:  mdh base (link 0) -> mdh frame (7)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S7RRRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 19:14:13
% EndTime: 2018-11-26 19:14:13
% DurationCPUTime: 0.16s
% Computational Cost: add. (242->52), mult. (600->71), div. (0->0), fcn. (857->14), ass. (0->50)
t37 = sin(qJ(3));
t38 = sin(qJ(2));
t52 = t38 * t37;
t44 = cos(qJ(3));
t51 = t38 * t44;
t39 = sin(qJ(1));
t50 = t39 * t38;
t45 = cos(qJ(2));
t49 = t39 * t45;
t46 = cos(qJ(1));
t48 = t46 * t38;
t47 = t46 * t45;
t32 = pkin(1) + 0;
t30 = t45 * pkin(2) + t32;
t26 = -pkin(2) * t50 + 0;
t27 = -pkin(2) * t48 + 0;
t36 = sin(qJ(4));
t43 = cos(qJ(4));
t20 = t36 * t51 + t45 * t43;
t19 = t20 * pkin(3) + t30;
t23 = t46 * t37 + t44 * t49;
t15 = t23 * t36 - t43 * t50;
t11 = t15 * pkin(3) + t26;
t25 = -t39 * t37 + t44 * t47;
t17 = t25 * t36 - t43 * t48;
t12 = t17 * pkin(3) + t27;
t42 = cos(qJ(5));
t41 = cos(qJ(6));
t40 = cos(qJ(7));
t35 = sin(qJ(5));
t34 = sin(qJ(6));
t33 = sin(qJ(7));
t24 = -t37 * t47 - t39 * t44;
t22 = -t37 * t49 + t46 * t44;
t21 = -t45 * t36 + t43 * t51;
t18 = t25 * t43 + t36 * t48;
t16 = t23 * t43 + t36 * t50;
t14 = t21 * t42 - t35 * t52;
t13 = t21 * t35 + t42 * t52;
t10 = t18 * t42 + t24 * t35;
t9 = t18 * t35 - t24 * t42;
t8 = t16 * t42 + t22 * t35;
t7 = t16 * t35 - t22 * t42;
t6 = t14 * t41 + t20 * t34;
t5 = -t14 * t34 + t20 * t41;
t4 = t10 * t41 + t17 * t34;
t3 = -t10 * t34 + t17 * t41;
t2 = t15 * t34 + t8 * t41;
t1 = t15 * t41 - t8 * t34;
t28 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t46, -t39, 0, 0; t39, t46, 0, 0; 0, 0, 1, t32; 0, 0, 0, 1; t47, -t48, t39, 0; t49, -t50, -t46, 0; t38, t45, 0, t32; 0, 0, 0, 1; t25, t24, -t48, t27; t23, t22, -t50, t26; t51, -t52, t45, t30; 0, 0, 0, 1; t18, -t17, t24, t27; t16, -t15, t22, t26; t21, -t20, -t52, t30; 0, 0, 0, 1; t10, -t9, t17, t12; t8, -t7, t15, t11; t14, -t13, t20, t19; 0, 0, 0, 1; t4, t3, t9, t12; t2, t1, t7, t11; t6, t5, t13, t19; 0, 0, 0, 1; -t9 * t33 + t4 * t40, -t4 * t33 - t9 * t40, t3, t3 * pkin(4) + t12; t2 * t40 - t7 * t33, -t2 * t33 - t7 * t40, t1, t1 * pkin(4) + t11; -t13 * t33 + t6 * t40, -t13 * t40 - t6 * t33, t5, t5 * pkin(4) + t19; 0, 0, 0, 1;];
T_ges = t28;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,7+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,7+1]); end % symbolisch
for i = 1:7+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
