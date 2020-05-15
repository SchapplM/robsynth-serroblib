% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2018-11-23 17:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPPR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:34:31
% EndTime: 2018-11-23 17:34:32
% DurationCPUTime: 0.16s
% Computational Cost: add. (204->64), mult. (230->71), div. (0->0), fcn. (314->10), ass. (0->39)
t24 = qJ(3) + pkin(10);
t19 = sin(t24);
t20 = cos(t24);
t34 = cos(qJ(1));
t30 = sin(qJ(1));
t33 = cos(qJ(2));
t48 = t30 * t33;
t3 = t19 * t48 + t34 * t20;
t4 = -t34 * t19 + t20 * t48;
t53 = t4 * pkin(4) + t3 * qJ(5);
t47 = t34 * t33;
t5 = t19 * t47 - t30 * t20;
t6 = t30 * t19 + t20 * t47;
t52 = t6 * pkin(4) + t5 * qJ(5);
t29 = sin(qJ(2));
t51 = t29 * t19;
t12 = t29 * t20;
t28 = sin(qJ(3));
t49 = t30 * t28;
t16 = t30 * t29;
t17 = t34 * t29;
t25 = pkin(6) + 0;
t45 = t30 * pkin(1) + 0;
t26 = -qJ(4) - pkin(8);
t44 = t29 * (-pkin(9) - t26);
t43 = t34 * pkin(1) + t30 * pkin(7) + 0;
t32 = cos(qJ(3));
t18 = t32 * pkin(3) + pkin(2);
t42 = t29 * t18 + t33 * t26 + t25;
t41 = pkin(2) * t33 + pkin(8) * t29;
t40 = -t34 * pkin(7) + t45;
t39 = pkin(3) * t49 + t18 * t47 + t43;
t38 = pkin(4) * t12 + qJ(5) * t51 + t42;
t37 = t18 * t48 + (-pkin(3) * t28 - pkin(7)) * t34 + t45;
t36 = -t26 * t17 + t39;
t35 = -t26 * t16 + t37;
t31 = cos(qJ(6));
t27 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t30, 0, 0; t30, t34, 0, 0; 0, 0, 1, t25; 0, 0, 0, 1; t47, -t17, t30, t43; t48, -t16, -t34, t40; t29, t33, 0, t25; 0, 0, 0, 1; t32 * t47 + t49, -t28 * t47 + t30 * t32, t17, t41 * t34 + t43; -t34 * t28 + t32 * t48, -t28 * t48 - t34 * t32, t16, t41 * t30 + t40; t29 * t32, -t29 * t28, -t33, t29 * pkin(2) - t33 * pkin(8) + t25; 0, 0, 0, 1; t6, -t5, t17, t36; t4, -t3, t16, t35; t12, -t51, -t33, t42; 0, 0, 0, 1; t6, t17, t5, t36 + t52; t4, t16, t3, t35 + t53; t12, -t33, t51, t38; 0, 0, 0, 1; t5 * t27 + t6 * t31, -t6 * t27 + t5 * t31, -t17, t6 * pkin(5) + t34 * t44 + t39 + t52; t3 * t27 + t4 * t31, -t4 * t27 + t3 * t31, -t16, t4 * pkin(5) + t30 * t44 + t37 + t53; (t19 * t27 + t20 * t31) * t29 (t19 * t31 - t20 * t27) * t29, t33, pkin(5) * t12 + t33 * pkin(9) + t38; 0, 0, 0, 1;];
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
