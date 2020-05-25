% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR2
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
% Datum: 2018-11-23 17:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:33:27
% EndTime: 2018-11-23 17:33:27
% DurationCPUTime: 0.14s
% Computational Cost: add. (186->55), mult. (100->50), div. (0->0), fcn. (149->10), ass. (0->35)
t29 = cos(qJ(1));
t23 = qJ(2) + qJ(3);
t16 = pkin(10) + t23;
t12 = sin(t16);
t38 = qJ(5) * t12;
t13 = cos(t16);
t9 = t29 * t13;
t45 = pkin(4) * t9 + t29 * t38;
t30 = -pkin(8) - pkin(7);
t28 = cos(qJ(2));
t15 = t28 * pkin(2) + pkin(1);
t26 = sin(qJ(1));
t44 = t26 * t12;
t8 = t26 * t13;
t24 = sin(qJ(6));
t43 = t26 * t24;
t27 = cos(qJ(6));
t42 = t26 * t27;
t41 = t29 * t12;
t40 = t29 * t24;
t39 = t29 * t27;
t22 = pkin(6) + 0;
t18 = cos(t23);
t3 = pkin(3) * t18 + t15;
t37 = t29 * t3 + 0;
t21 = -qJ(4) + t30;
t36 = t29 * t21 + t26 * t3 + 0;
t25 = sin(qJ(2));
t35 = t25 * pkin(2) + t22;
t17 = sin(t23);
t34 = pkin(3) * t17 + t35;
t33 = pkin(4) * t8 + t26 * t38 + t36;
t32 = -t26 * t21 + t37;
t31 = t12 * pkin(4) - t13 * qJ(5) + t34;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t29 * t28, -t29 * t25, t26, t29 * pkin(1) + t26 * pkin(7) + 0; t26 * t28, -t26 * t25, -t29, t26 * pkin(1) - t29 * pkin(7) + 0; t25, t28, 0, t22; 0, 0, 0, 1; t29 * t18, -t29 * t17, t26, t29 * t15 - t26 * t30 + 0; t26 * t18, -t26 * t17, -t29, t26 * t15 + t29 * t30 + 0; t17, t18, 0, t35; 0, 0, 0, 1; t9, -t41, t26, t32; t8, -t44, -t29, t36; t12, t13, 0, t34; 0, 0, 0, 1; t26, -t9, t41, t32 + t45; -t29, -t8, t44, t33; 0, -t12, -t13, t31; 0, 0, 0, 1; t12 * t40 + t42, t12 * t39 - t43, t9, pkin(9) * t9 + (pkin(5) - t21) * t26 + t37 + t45; t12 * t43 - t39, t12 * t42 + t40, t8, -t29 * pkin(5) + pkin(9) * t8 + t33; -t13 * t24, -t13 * t27, t12, t12 * pkin(9) + t31; 0, 0, 0, 1;];
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
