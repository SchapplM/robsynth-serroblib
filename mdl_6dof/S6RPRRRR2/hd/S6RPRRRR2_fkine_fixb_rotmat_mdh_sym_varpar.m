% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:32:45
% EndTime: 2018-11-23 16:32:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (200->53), mult. (102->60), div. (0->0), fcn. (152->12), ass. (0->37)
t21 = sin(qJ(5));
t18 = qJ(1) + pkin(11);
t9 = cos(t18);
t42 = t9 * t21;
t19 = qJ(5) + qJ(6);
t11 = sin(t19);
t20 = qJ(3) + qJ(4);
t14 = cos(t20);
t41 = t11 * t14;
t13 = cos(t19);
t40 = t13 * t14;
t39 = t14 * t21;
t24 = cos(qJ(5));
t38 = t14 * t24;
t37 = pkin(6) + 0;
t23 = sin(qJ(1));
t36 = t23 * pkin(1) + 0;
t26 = cos(qJ(1));
t35 = t26 * pkin(1) + 0;
t25 = cos(qJ(3));
t7 = t25 * pkin(3) + pkin(2);
t34 = t9 * t7 + t35;
t10 = qJ(2) + t37;
t28 = -pkin(8) - pkin(7);
t8 = sin(t18);
t33 = t9 * t28 + t8 * t7 + t36;
t22 = sin(qJ(3));
t32 = t22 * pkin(3) + t10;
t12 = sin(t20);
t31 = pkin(4) * t14 + pkin(9) * t12;
t27 = -pkin(10) - pkin(9);
t6 = t24 * pkin(5) + pkin(4);
t30 = -t12 * t27 + t14 * t6;
t29 = -t8 * t28 + t34;
t4 = t9 * t12;
t3 = t8 * t12;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t23, 0, 0; t23, t26, 0, 0; 0, 0, 1, t37; 0, 0, 0, 1; t9, -t8, 0, t35; t8, t9, 0, t36; 0, 0, 1, t10; 0, 0, 0, 1; t9 * t25, -t9 * t22, t8, t9 * pkin(2) + t8 * pkin(7) + t35; t8 * t25, -t8 * t22, -t9, t8 * pkin(2) - t9 * pkin(7) + t36; t22, t25, 0, t10; 0, 0, 0, 1; t9 * t14, -t4, t8, t29; t8 * t14, -t3, -t9, t33; t12, t14, 0, t32; 0, 0, 0, 1; t8 * t21 + t9 * t38, t8 * t24 - t9 * t39, t4, t31 * t9 + t29; t8 * t38 - t42, -t9 * t24 - t8 * t39, t3, t31 * t8 + t33; t12 * t24, -t12 * t21, -t14, t12 * pkin(4) - t14 * pkin(9) + t32; 0, 0, 0, 1; t8 * t11 + t9 * t40, t8 * t13 - t9 * t41, t4, t30 * t9 + (pkin(5) * t21 - t28) * t8 + t34; -t9 * t11 + t8 * t40, -t9 * t13 - t8 * t41, t3, -pkin(5) * t42 + t30 * t8 + t33; t12 * t13, -t12 * t11, -t14, t12 * t6 + t14 * t27 + t32; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
