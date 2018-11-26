% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRP12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:19:04
% EndTime: 2018-11-23 17:19:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (147->62), mult. (177->57), div. (0->0), fcn. (240->8), ass. (0->40)
t26 = sin(qJ(1));
t28 = cos(qJ(2));
t13 = t26 * t28;
t25 = sin(qJ(2));
t44 = qJ(3) * t25;
t51 = pkin(2) * t13 + t26 * t44;
t50 = t26 * t25;
t27 = cos(qJ(4));
t49 = t26 * t27;
t23 = qJ(4) + qJ(5);
t16 = sin(t23);
t48 = t28 * t16;
t17 = cos(t23);
t47 = t28 * t17;
t29 = cos(qJ(1));
t46 = t29 * t25;
t45 = t29 * t27;
t14 = t29 * t28;
t22 = pkin(6) + 0;
t43 = t26 * pkin(1) + 0;
t24 = sin(qJ(4));
t42 = t24 * t50;
t41 = t24 * t46;
t40 = t25 * pkin(2) + t22;
t39 = -pkin(4) * t24 - qJ(3);
t38 = t29 * pkin(1) + t26 * pkin(7) + 0;
t37 = t43 + t51;
t36 = -t29 * pkin(7) + t43;
t35 = pkin(2) * t14 + t29 * t44 + t38;
t30 = -pkin(9) - pkin(8);
t34 = -t25 * t30 + t40;
t33 = -t28 * qJ(3) + t40;
t15 = t27 * pkin(4) + pkin(3);
t32 = pkin(4) * t41 - t30 * t14 + t26 * t15 + t35;
t31 = -t30 * t13 + pkin(4) * t42 + (-pkin(7) - t15) * t29 + t37;
t4 = t16 * t50 - t29 * t17;
t3 = t29 * t16 + t17 * t50;
t2 = t16 * t46 + t26 * t17;
t1 = t26 * t16 - t17 * t46;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t14, -t46, t26, t38; t13, -t50, -t29, t36; t25, t28, 0, t22; 0, 0, 0, 1; t26, -t14, t46, t35; -t29, -t13, t50, t36 + t51; 0, -t25, -t28, t33; 0, 0, 0, 1; t41 + t49, -t26 * t24 + t25 * t45, t14, t26 * pkin(3) + pkin(8) * t14 + t35; t42 - t45, t29 * t24 + t25 * t49, t13, pkin(8) * t13 + (-pkin(3) - pkin(7)) * t29 + t37; -t28 * t24, -t28 * t27, t25, t25 * pkin(8) + t33; 0, 0, 0, 1; t2, -t1, t14, t32; t4, t3, t13, t31; -t48, -t47, t25, t39 * t28 + t34; 0, 0, 0, 1; t2, t14, t1, t2 * pkin(5) + t1 * qJ(6) + t32; t4, t13, -t3, t4 * pkin(5) - t3 * qJ(6) + t31; -t48, t25, t47 (-pkin(5) * t16 + qJ(6) * t17 + t39) * t28 + t34; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
