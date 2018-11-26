% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRPR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:13:47
% EndTime: 2018-11-23 18:13:47
% DurationCPUTime: 0.15s
% Computational Cost: add. (199->66), mult. (137->74), div. (0->0), fcn. (196->12), ass. (0->40)
t24 = sin(qJ(4));
t44 = t24 * pkin(4);
t27 = cos(qJ(4));
t9 = t27 * pkin(4) + pkin(3);
t22 = qJ(2) + qJ(3);
t16 = cos(t22);
t26 = sin(qJ(1));
t43 = t26 * t16;
t42 = t26 * t24;
t41 = t26 * t27;
t29 = cos(qJ(1));
t40 = t29 * t16;
t39 = t29 * t24;
t38 = t29 * t27;
t23 = -qJ(5) - pkin(9);
t21 = pkin(6) + 0;
t20 = qJ(4) + pkin(11);
t28 = cos(qJ(2));
t10 = t28 * pkin(2) + pkin(1);
t37 = t29 * t10 + 0;
t25 = sin(qJ(2));
t36 = t25 * pkin(2) + t21;
t30 = -pkin(8) - pkin(7);
t35 = t26 * t10 + t29 * t30 + 0;
t15 = sin(t22);
t34 = pkin(3) * t16 + pkin(9) * t15;
t13 = cos(t20);
t1 = pkin(5) * t13 + t9;
t19 = -pkin(10) + t23;
t33 = t1 * t16 - t15 * t19;
t32 = -t15 * t23 + t16 * t9;
t31 = -t26 * t30 + t37;
t14 = qJ(6) + t20;
t12 = sin(t20);
t8 = cos(t14);
t7 = sin(t14);
t6 = t29 * t15;
t5 = t26 * t15;
t2 = pkin(5) * t12 + t44;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t29 * t28, -t29 * t25, t26, t29 * pkin(1) + t26 * pkin(7) + 0; t26 * t28, -t26 * t25, -t29, t26 * pkin(1) - t29 * pkin(7) + 0; t25, t28, 0, t21; 0, 0, 0, 1; t40, -t6, t26, t31; t43, -t5, -t29, t35; t15, t16, 0, t36; 0, 0, 0, 1; t16 * t38 + t42, -t16 * t39 + t41, t6, t34 * t29 + t31; t16 * t41 - t39, -t16 * t42 - t38, t5, t34 * t26 + t35; t15 * t27, -t15 * t24, -t16, t15 * pkin(3) - t16 * pkin(9) + t36; 0, 0, 0, 1; t26 * t12 + t13 * t40, -t12 * t40 + t26 * t13, t6, t32 * t29 + (-t30 + t44) * t26 + t37; -t29 * t12 + t13 * t43, -t12 * t43 - t29 * t13, t5, -pkin(4) * t39 + t32 * t26 + t35; t15 * t13, -t15 * t12, -t16, t15 * t9 + t16 * t23 + t36; 0, 0, 0, 1; t26 * t7 + t8 * t40, t26 * t8 - t7 * t40, t6, t33 * t29 + (t2 - t30) * t26 + t37; -t29 * t7 + t8 * t43, -t29 * t8 - t7 * t43, t5, -t29 * t2 + t33 * t26 + t35; t15 * t8, -t15 * t7, -t16, t15 * t1 + t16 * t19 + t36; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
