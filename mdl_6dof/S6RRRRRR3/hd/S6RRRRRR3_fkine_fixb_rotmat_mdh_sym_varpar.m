% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:38:47
% EndTime: 2018-11-23 18:38:48
% DurationCPUTime: 0.15s
% Computational Cost: add. (199->66), mult. (137->74), div. (0->0), fcn. (196->12), ass. (0->40)
t29 = -pkin(10) - pkin(9);
t23 = sin(qJ(4));
t44 = t23 * pkin(4);
t26 = cos(qJ(4));
t9 = t26 * pkin(4) + pkin(3);
t22 = qJ(2) + qJ(3);
t15 = cos(t22);
t25 = sin(qJ(1));
t43 = t25 * t15;
t42 = t25 * t23;
t41 = t25 * t26;
t28 = cos(qJ(1));
t40 = t28 * t15;
t39 = t28 * t23;
t38 = t28 * t26;
t21 = qJ(4) + qJ(5);
t19 = pkin(6) + 0;
t27 = cos(qJ(2));
t10 = t27 * pkin(2) + pkin(1);
t37 = t28 * t10 + 0;
t24 = sin(qJ(2));
t36 = t24 * pkin(2) + t19;
t30 = -pkin(8) - pkin(7);
t35 = t25 * t10 + t28 * t30 + 0;
t13 = sin(t22);
t34 = pkin(3) * t15 + pkin(9) * t13;
t14 = cos(t21);
t1 = pkin(5) * t14 + t9;
t20 = -pkin(11) + t29;
t33 = t1 * t15 - t13 * t20;
t32 = -t13 * t29 + t15 * t9;
t31 = -t25 * t30 + t37;
t16 = qJ(6) + t21;
t12 = sin(t21);
t8 = cos(t16);
t7 = sin(t16);
t6 = t28 * t13;
t5 = t25 * t13;
t2 = pkin(5) * t12 + t44;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t19; 0, 0, 0, 1; t28 * t27, -t28 * t24, t25, t28 * pkin(1) + t25 * pkin(7) + 0; t25 * t27, -t25 * t24, -t28, t25 * pkin(1) - t28 * pkin(7) + 0; t24, t27, 0, t19; 0, 0, 0, 1; t40, -t6, t25, t31; t43, -t5, -t28, t35; t13, t15, 0, t36; 0, 0, 0, 1; t15 * t38 + t42, -t15 * t39 + t41, t6, t34 * t28 + t31; t15 * t41 - t39, -t15 * t42 - t38, t5, t34 * t25 + t35; t13 * t26, -t13 * t23, -t15, t13 * pkin(3) - t15 * pkin(9) + t36; 0, 0, 0, 1; t25 * t12 + t14 * t40, -t12 * t40 + t25 * t14, t6, t32 * t28 + (-t30 + t44) * t25 + t37; -t28 * t12 + t14 * t43, -t12 * t43 - t28 * t14, t5, -pkin(4) * t39 + t32 * t25 + t35; t13 * t14, -t13 * t12, -t15, t13 * t9 + t15 * t29 + t36; 0, 0, 0, 1; t25 * t7 + t8 * t40, t25 * t8 - t7 * t40, t6, t33 * t28 + (t2 - t30) * t25 + t37; -t28 * t7 + t8 * t43, -t28 * t8 - t7 * t43, t5, -t28 * t2 + t33 * t25 + t35; t13 * t8, -t13 * t7, -t15, t13 * t1 + t15 * t20 + t36; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
