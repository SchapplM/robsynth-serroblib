% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2018-11-23 16:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRPP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:11:10
% EndTime: 2018-11-23 16:11:10
% DurationCPUTime: 0.12s
% Computational Cost: add. (213->52), mult. (142->56), div. (0->0), fcn. (201->10), ass. (0->40)
t26 = qJ(1) + pkin(9);
t19 = sin(t26);
t28 = sin(qJ(4));
t50 = t19 * t28;
t32 = cos(qJ(3));
t49 = t19 * t32;
t21 = cos(t26);
t48 = t21 * t32;
t27 = -qJ(5) - pkin(8);
t29 = sin(qJ(3));
t47 = t27 * t29;
t46 = t28 * t32;
t25 = qJ(4) + pkin(10);
t18 = sin(t25);
t45 = t29 * t18;
t31 = cos(qJ(4));
t44 = t31 * t32;
t43 = pkin(6) + 0;
t30 = sin(qJ(1));
t42 = t30 * pkin(1) + 0;
t33 = cos(qJ(1));
t41 = t33 * pkin(1) + 0;
t22 = qJ(2) + t43;
t40 = t19 * pkin(2) + t42;
t39 = t21 * pkin(2) + t19 * pkin(7) + t41;
t38 = pkin(3) * t32 + pkin(8) * t29;
t17 = t31 * pkin(4) + pkin(3);
t37 = t29 * t17 + t32 * t27 + t22;
t36 = -t21 * pkin(7) + t40;
t35 = pkin(4) * t50 + t17 * t48 - t21 * t47 + t39;
t34 = -t19 * t47 + t17 * t49 + (-pkin(4) * t28 - pkin(7)) * t21 + t40;
t20 = cos(t25);
t12 = t21 * t29;
t11 = t29 * t20;
t10 = t19 * t29;
t4 = t19 * t18 + t20 * t48;
t3 = t18 * t48 - t19 * t20;
t2 = -t21 * t18 + t20 * t49;
t1 = t18 * t49 + t21 * t20;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t30, 0, 0; t30, t33, 0, 0; 0, 0, 1, t43; 0, 0, 0, 1; t21, -t19, 0, t41; t19, t21, 0, t42; 0, 0, 1, t22; 0, 0, 0, 1; t48, -t12, t19, t39; t49, -t10, -t21, t36; t29, t32, 0, t22; 0, 0, 0, 1; t21 * t44 + t50, t19 * t31 - t21 * t46, t12, t38 * t21 + t39; t19 * t44 - t21 * t28, -t19 * t46 - t21 * t31, t10, t38 * t19 + t36; t29 * t31, -t29 * t28, -t32, t29 * pkin(3) - t32 * pkin(8) + t22; 0, 0, 0, 1; t4, -t3, t12, t35; t2, -t1, t10, t34; t11, -t45, -t32, t37; 0, 0, 0, 1; t4, t12, t3, t4 * pkin(5) + t3 * qJ(6) + t35; t2, t10, t1, t2 * pkin(5) + t1 * qJ(6) + t34; t11, -t32, t45 (pkin(5) * t20 + qJ(6) * t18) * t29 + t37; 0, 0, 0, 1;];
T_ges = t5;
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
