% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR1
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
% Datum: 2018-11-23 16:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:02:54
% EndTime: 2018-11-23 16:02:54
% DurationCPUTime: 0.13s
% Computational Cost: add. (191->47), mult. (79->46), div. (0->0), fcn. (124->12), ass. (0->35)
t21 = qJ(1) + pkin(10);
t10 = sin(t21);
t20 = qJ(3) + pkin(11);
t13 = qJ(5) + t20;
t6 = sin(t13);
t42 = t10 * t6;
t12 = cos(t21);
t41 = t12 * t6;
t27 = cos(qJ(3));
t8 = t27 * pkin(3) + pkin(2);
t23 = sin(qJ(6));
t40 = t10 * t23;
t26 = cos(qJ(6));
t39 = t10 * t26;
t38 = t12 * t23;
t37 = t12 * t26;
t22 = -qJ(4) - pkin(7);
t36 = pkin(6) + 0;
t25 = sin(qJ(1));
t35 = t25 * pkin(1) + 0;
t28 = cos(qJ(1));
t34 = t28 * pkin(1) + 0;
t14 = qJ(2) + t36;
t19 = -pkin(8) + t22;
t11 = cos(t20);
t3 = pkin(4) * t11 + t8;
t33 = t10 * t3 + t12 * t19 + t35;
t7 = cos(t13);
t32 = pkin(5) * t7 + pkin(9) * t6;
t24 = sin(qJ(3));
t31 = t24 * pkin(3) + t14;
t9 = sin(t20);
t30 = pkin(4) * t9 + t31;
t29 = -t10 * t19 + t12 * t3 + t34;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t36; 0, 0, 0, 1; t12, -t10, 0, t34; t10, t12, 0, t35; 0, 0, 1, t14; 0, 0, 0, 1; t12 * t27, -t12 * t24, t10, t12 * pkin(2) + t10 * pkin(7) + t34; t10 * t27, -t10 * t24, -t12, t10 * pkin(2) - t12 * pkin(7) + t35; t24, t27, 0, t14; 0, 0, 0, 1; t12 * t11, -t12 * t9, t10, -t10 * t22 + t12 * t8 + t34; t10 * t11, -t10 * t9, -t12, t10 * t8 + t12 * t22 + t35; t9, t11, 0, t31; 0, 0, 0, 1; t12 * t7, -t41, t10, t29; t10 * t7, -t42, -t12, t33; t6, t7, 0, t30; 0, 0, 0, 1; t7 * t37 + t40, -t7 * t38 + t39, t41, t32 * t12 + t29; t7 * t39 - t38, -t7 * t40 - t37, t42, t32 * t10 + t33; t6 * t26, -t6 * t23, -t7, t6 * pkin(5) - t7 * pkin(9) + t30; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
