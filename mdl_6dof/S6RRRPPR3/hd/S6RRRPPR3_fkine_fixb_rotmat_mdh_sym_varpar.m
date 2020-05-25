% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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

function [T_c_mdh, Tc_stack] = S6RRRPPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:33:56
% EndTime: 2018-11-23 17:33:56
% DurationCPUTime: 0.11s
% Computational Cost: add. (156->52), mult. (118->46), div. (0->0), fcn. (167->8), ass. (0->34)
t22 = qJ(2) + qJ(3);
t18 = cos(t22);
t28 = cos(qJ(1));
t12 = t28 * t18;
t17 = sin(t22);
t40 = qJ(4) * t17;
t45 = pkin(3) * t12 + t28 * t40;
t25 = sin(qJ(1));
t10 = t25 * t18;
t23 = sin(qJ(6));
t44 = t25 * t23;
t26 = cos(qJ(6));
t43 = t25 * t26;
t42 = t28 * t23;
t41 = t28 * t26;
t21 = pkin(6) + 0;
t27 = cos(qJ(2));
t15 = t27 * pkin(2) + pkin(1);
t39 = t28 * t15 + 0;
t24 = sin(qJ(2));
t38 = t24 * pkin(2) + t21;
t29 = -pkin(8) - pkin(7);
t37 = t25 * t15 + t28 * t29 + 0;
t36 = t17 * pkin(3) + t38;
t35 = pkin(5) * t17 + pkin(9) * t18;
t34 = pkin(3) * t10 + t25 * t40 + t37;
t33 = -t25 * t29 + t39;
t32 = pkin(4) * t10 + t28 * qJ(5) + t34;
t31 = -t18 * qJ(4) + t36;
t30 = pkin(4) * t12 + (-qJ(5) - t29) * t25 + t39 + t45;
t13 = t17 * pkin(4);
t11 = t28 * t17;
t9 = t25 * t17;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t28 * t27, -t28 * t24, t25, t28 * pkin(1) + t25 * pkin(7) + 0; t25 * t27, -t25 * t24, -t28, t25 * pkin(1) - t28 * pkin(7) + 0; t24, t27, 0, t21; 0, 0, 0, 1; t12, -t11, t25, t33; t10, -t9, -t28, t37; t17, t18, 0, t38; 0, 0, 0, 1; t12, t25, t11, t33 + t45; t10, -t28, t9, t34; t17, 0, -t18, t31; 0, 0, 0, 1; t11, -t12, -t25, t30; t9, -t10, t28, t32; -t18, -t17, 0, t13 + t31; 0, 0, 0, 1; t17 * t41 - t44, -t17 * t42 - t43, t12, t35 * t28 + t30; t17 * t43 + t42, -t17 * t44 + t41, t10, t35 * t25 + t32; -t18 * t26, t18 * t23, t17, t17 * pkin(9) + t13 + (-pkin(5) - qJ(4)) * t18 + t36; 0, 0, 0, 1;];
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
