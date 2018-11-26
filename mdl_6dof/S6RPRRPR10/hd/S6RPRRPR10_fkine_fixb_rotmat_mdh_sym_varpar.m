% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2018-11-23 16:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:21:21
% EndTime: 2018-11-23 16:21:21
% DurationCPUTime: 0.14s
% Computational Cost: add. (118->56), mult. (191->54), div. (0->0), fcn. (266->8), ass. (0->37)
t24 = sin(qJ(4));
t25 = sin(qJ(3));
t30 = cos(qJ(1));
t46 = t30 * t25;
t26 = sin(qJ(1));
t28 = cos(qJ(4));
t48 = t26 * t28;
t5 = t24 * t46 + t48;
t45 = t30 * t28;
t6 = t26 * t24 - t25 * t45;
t50 = t6 * pkin(4) - t5 * qJ(5);
t49 = t26 * t25;
t29 = cos(qJ(3));
t11 = t26 * t29;
t47 = t29 * t24;
t12 = t29 * t28;
t14 = t30 * t29;
t22 = pkin(6) + 0;
t43 = pkin(8) * t11;
t42 = t26 * pkin(1) + 0;
t41 = pkin(2) + t22;
t40 = -pkin(3) * t25 - qJ(2);
t39 = t30 * pkin(1) + t26 * qJ(2) + 0;
t17 = t26 * pkin(7);
t38 = pkin(8) * t14 + t17 + t42;
t37 = t30 * pkin(7) + t39;
t36 = t29 * pkin(3) + t25 * pkin(8) + t41;
t35 = pkin(3) * t49 + t37;
t34 = -t30 * qJ(2) + t42;
t33 = pkin(4) * t12 + qJ(5) * t47 + t36;
t3 = t24 * t49 - t45;
t4 = t30 * t24 + t25 * t48;
t32 = t4 * pkin(4) + t3 * qJ(5) + t35;
t31 = t40 * t30 + t38;
t27 = cos(qJ(6));
t23 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t26, 0, 0; t26, t30, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; 0, -t30, t26, t39; 0, -t26, -t30, t34; 1, 0, 0, t22; 0, 0, 0, 1; t49, t11, t30, t37; -t46, -t14, t26, t17 + t34; t29, -t25, 0, t41; 0, 0, 0, 1; t4, -t3, -t11, t35 - t43; t6, t5, t14, t31; t12, -t47, t25, t36; 0, 0, 0, 1; t4, -t11, t3, t32 - t43; t6, t14, -t5, t31 + t50; t12, t25, t47, t33; 0, 0, 0, 1; t3 * t23 + t4 * t27, -t4 * t23 + t3 * t27, t11, t4 * pkin(5) + (-pkin(8) + pkin(9)) * t11 + t32; -t5 * t23 + t6 * t27, -t6 * t23 - t5 * t27, -t14, t6 * pkin(5) + (-pkin(9) * t29 + t40) * t30 + t38 + t50; (t23 * t24 + t27 * t28) * t29 (-t23 * t28 + t24 * t27) * t29, -t25, pkin(5) * t12 - t25 * pkin(9) + t33; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
