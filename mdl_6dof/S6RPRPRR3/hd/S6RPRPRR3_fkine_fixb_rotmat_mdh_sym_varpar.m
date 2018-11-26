% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2018-11-23 16:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRPRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:03:47
% EndTime: 2018-11-23 16:03:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (209->60), mult. (127->68), div. (0->0), fcn. (182->12), ass. (0->39)
t23 = sin(pkin(11));
t44 = t23 * pkin(4);
t24 = cos(pkin(11));
t10 = t24 * pkin(4) + pkin(3);
t22 = qJ(1) + pkin(10);
t12 = sin(t22);
t43 = t12 * t23;
t28 = cos(qJ(3));
t42 = t12 * t28;
t14 = cos(t22);
t41 = t14 * t28;
t40 = t23 * t28;
t39 = t24 * t28;
t25 = -pkin(8) - qJ(4);
t38 = pkin(6) + 0;
t21 = pkin(11) + qJ(5);
t27 = sin(qJ(1));
t37 = t27 * pkin(1) + 0;
t29 = cos(qJ(1));
t36 = t29 * pkin(1) + 0;
t35 = t12 * pkin(2) + t37;
t16 = qJ(2) + t38;
t34 = t14 * pkin(2) + t12 * pkin(7) + t36;
t13 = cos(t21);
t1 = pkin(5) * t13 + t10;
t20 = -pkin(9) + t25;
t26 = sin(qJ(3));
t33 = t1 * t28 - t20 * t26;
t32 = pkin(3) * t28 + qJ(4) * t26;
t31 = t10 * t28 - t25 * t26;
t30 = -t14 * pkin(7) + t35;
t15 = qJ(6) + t21;
t11 = sin(t21);
t9 = cos(t15);
t8 = sin(t15);
t4 = t14 * t26;
t3 = t12 * t26;
t2 = pkin(5) * t11 + t44;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t27, 0, 0; t27, t29, 0, 0; 0, 0, 1, t38; 0, 0, 0, 1; t14, -t12, 0, t36; t12, t14, 0, t37; 0, 0, 1, t16; 0, 0, 0, 1; t41, -t4, t12, t34; t42, -t3, -t14, t30; t26, t28, 0, t16; 0, 0, 0, 1; t14 * t39 + t43, t12 * t24 - t14 * t40, t4, t32 * t14 + t34; t12 * t39 - t14 * t23, -t12 * t40 - t14 * t24, t3, t32 * t12 + t30; t26 * t24, -t26 * t23, -t28, t26 * pkin(3) - t28 * qJ(4) + t16; 0, 0, 0, 1; t12 * t11 + t13 * t41, -t11 * t41 + t12 * t13, t4, pkin(4) * t43 + t31 * t14 + t34; -t14 * t11 + t13 * t42, -t11 * t42 - t14 * t13, t3 (-pkin(7) - t44) * t14 + t31 * t12 + t35; t26 * t13, -t26 * t11, -t28, t26 * t10 + t28 * t25 + t16; 0, 0, 0, 1; t12 * t8 + t9 * t41, t12 * t9 - t8 * t41, t4, t12 * t2 + t33 * t14 + t34; -t14 * t8 + t9 * t42, -t14 * t9 - t8 * t42, t3 (-pkin(7) - t2) * t14 + t33 * t12 + t35; t26 * t9, -t26 * t8, -t28, t26 * t1 + t28 * t20 + t16; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
