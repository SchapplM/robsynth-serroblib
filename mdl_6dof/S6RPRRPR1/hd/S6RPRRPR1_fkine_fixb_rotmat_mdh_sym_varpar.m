% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 16:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:15:41
% EndTime: 2018-11-23 16:15:41
% DurationCPUTime: 0.13s
% Computational Cost: add. (191->47), mult. (79->46), div. (0->0), fcn. (124->12), ass. (0->35)
t21 = qJ(3) + qJ(4);
t11 = pkin(11) + t21;
t5 = sin(t11);
t20 = qJ(1) + pkin(10);
t9 = sin(t20);
t42 = t9 * t5;
t28 = -pkin(8) - pkin(7);
t10 = cos(t20);
t41 = t10 * t5;
t22 = sin(qJ(6));
t40 = t9 * t22;
t25 = cos(qJ(6));
t39 = t9 * t25;
t26 = cos(qJ(3));
t8 = t26 * pkin(3) + pkin(2);
t38 = t10 * t22;
t37 = t10 * t25;
t36 = pkin(6) + 0;
t24 = sin(qJ(1));
t35 = t24 * pkin(1) + 0;
t27 = cos(qJ(1));
t34 = t27 * pkin(1) + 0;
t12 = qJ(2) + t36;
t19 = -qJ(5) + t28;
t14 = cos(t21);
t3 = pkin(4) * t14 + t8;
t33 = t10 * t19 + t9 * t3 + t35;
t6 = cos(t11);
t32 = pkin(5) * t6 + pkin(9) * t5;
t23 = sin(qJ(3));
t31 = t23 * pkin(3) + t12;
t13 = sin(t21);
t30 = pkin(4) * t13 + t31;
t29 = t10 * t3 - t9 * t19 + t34;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t24, 0, 0; t24, t27, 0, 0; 0, 0, 1, t36; 0, 0, 0, 1; t10, -t9, 0, t34; t9, t10, 0, t35; 0, 0, 1, t12; 0, 0, 0, 1; t10 * t26, -t10 * t23, t9, t10 * pkin(2) + t9 * pkin(7) + t34; t9 * t26, -t9 * t23, -t10, t9 * pkin(2) - t10 * pkin(7) + t35; t23, t26, 0, t12; 0, 0, 0, 1; t10 * t14, -t10 * t13, t9, t10 * t8 - t9 * t28 + t34; t9 * t14, -t9 * t13, -t10, t10 * t28 + t9 * t8 + t35; t13, t14, 0, t31; 0, 0, 0, 1; t10 * t6, -t41, t9, t29; t9 * t6, -t42, -t10, t33; t5, t6, 0, t30; 0, 0, 0, 1; t37 * t6 + t40, -t38 * t6 + t39, t41, t10 * t32 + t29; t39 * t6 - t38, -t40 * t6 - t37, t42, t32 * t9 + t33; t5 * t25, -t5 * t22, -t6, t5 * pkin(5) - t6 * pkin(9) + t30; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
