% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:43:22
% EndTime: 2018-11-23 17:43:22
% DurationCPUTime: 0.15s
% Computational Cost: add. (207->65), mult. (179->71), div. (0->0), fcn. (247->10), ass. (0->39)
t30 = sin(qJ(3));
t48 = t30 * pkin(3);
t33 = cos(qJ(3));
t18 = t33 * pkin(3) + pkin(2);
t27 = qJ(3) + pkin(10);
t21 = qJ(5) + t27;
t14 = sin(t21);
t31 = sin(qJ(2));
t47 = t31 * t14;
t32 = sin(qJ(1));
t46 = t32 * t30;
t16 = t32 * t31;
t34 = cos(qJ(2));
t45 = t32 * t34;
t35 = cos(qJ(1));
t17 = t35 * t31;
t44 = t35 * t34;
t29 = -qJ(4) - pkin(8);
t28 = pkin(6) + 0;
t43 = t32 * pkin(1) + 0;
t42 = t35 * pkin(1) + t32 * pkin(7) + 0;
t26 = -pkin(9) + t29;
t20 = cos(t27);
t9 = pkin(4) * t20 + t18;
t41 = t34 * t26 + t31 * t9 + t28;
t40 = pkin(2) * t34 + pkin(8) * t31;
t39 = t18 * t34 - t29 * t31;
t38 = -t35 * pkin(7) + t43;
t19 = sin(t27);
t10 = pkin(4) * t19 + t48;
t37 = t32 * t10 - t26 * t17 + t9 * t44 + t42;
t36 = -t26 * t16 + t9 * t45 + (-pkin(7) - t10) * t35 + t43;
t15 = cos(t21);
t11 = t31 * t15;
t4 = t32 * t14 + t15 * t44;
t3 = t14 * t44 - t32 * t15;
t2 = -t35 * t14 + t15 * t45;
t1 = t14 * t45 + t35 * t15;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t32, 0, 0; t32, t35, 0, 0; 0, 0, 1, t28; 0, 0, 0, 1; t44, -t17, t32, t42; t45, -t16, -t35, t38; t31, t34, 0, t28; 0, 0, 0, 1; t33 * t44 + t46, -t30 * t44 + t32 * t33, t17, t40 * t35 + t42; -t35 * t30 + t33 * t45, -t30 * t45 - t35 * t33, t16, t40 * t32 + t38; t31 * t33, -t31 * t30, -t34, t31 * pkin(2) - t34 * pkin(8) + t28; 0, 0, 0, 1; t32 * t19 + t20 * t44, -t19 * t44 + t32 * t20, t17, pkin(3) * t46 + t39 * t35 + t42; -t35 * t19 + t20 * t45, -t19 * t45 - t35 * t20, t16 (-pkin(7) - t48) * t35 + t39 * t32 + t43; t31 * t20, -t31 * t19, -t34, t31 * t18 + t34 * t29 + t28; 0, 0, 0, 1; t4, -t3, t17, t37; t2, -t1, t16, t36; t11, -t47, -t34, t41; 0, 0, 0, 1; t4, t17, t3, t4 * pkin(5) + t3 * qJ(6) + t37; t2, t16, t1, t2 * pkin(5) + t1 * qJ(6) + t36; t11, -t34, t47 (pkin(5) * t15 + qJ(6) * t14) * t31 + t41; 0, 0, 0, 1;];
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
