% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2018-11-23 15:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:59:15
% EndTime: 2018-11-23 15:59:15
% DurationCPUTime: 0.12s
% Computational Cost: add. (201->57), mult. (152->60), div. (0->0), fcn. (215->10), ass. (0->39)
t24 = pkin(10) + qJ(5);
t19 = sin(t24);
t25 = pkin(9) + qJ(3);
t20 = sin(t25);
t48 = t20 * t19;
t33 = sin(qJ(1));
t12 = t33 * t20;
t22 = cos(t25);
t47 = t33 * t22;
t27 = sin(pkin(10));
t46 = t33 * t27;
t29 = cos(pkin(10));
t45 = t33 * t29;
t34 = cos(qJ(1));
t14 = t34 * t20;
t44 = t34 * t22;
t43 = t34 * t27;
t42 = t34 * t29;
t26 = pkin(6) + 0;
t28 = sin(pkin(9));
t41 = t28 * pkin(2) + t26;
t30 = cos(pkin(9));
t17 = t30 * pkin(2) + pkin(1);
t32 = -pkin(7) - qJ(2);
t40 = t33 * t17 + t34 * t32 + 0;
t39 = pkin(3) * t22 + qJ(4) * t20;
t16 = t29 * pkin(4) + pkin(3);
t31 = -pkin(8) - qJ(4);
t38 = t20 * t16 + t22 * t31 + t41;
t37 = t34 * t17 - t33 * t32 + 0;
t36 = pkin(4) * t46 - t31 * t14 + t16 * t44 + t37;
t35 = -pkin(4) * t43 - t31 * t12 + t16 * t47 + t40;
t21 = cos(t24);
t8 = t20 * t21;
t4 = t33 * t19 + t21 * t44;
t3 = t19 * t44 - t33 * t21;
t2 = -t34 * t19 + t21 * t47;
t1 = t19 * t47 + t34 * t21;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t33, 0, 0; t33, t34, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t34 * t30, -t34 * t28, t33, t34 * pkin(1) + t33 * qJ(2) + 0; t33 * t30, -t33 * t28, -t34, t33 * pkin(1) - t34 * qJ(2) + 0; t28, t30, 0, t26; 0, 0, 0, 1; t44, -t14, t33, t37; t47, -t12, -t34, t40; t20, t22, 0, t41; 0, 0, 0, 1; t22 * t42 + t46, -t22 * t43 + t45, t14, t39 * t34 + t37; t22 * t45 - t43, -t22 * t46 - t42, t12, t39 * t33 + t40; t20 * t29, -t20 * t27, -t22, t20 * pkin(3) - t22 * qJ(4) + t41; 0, 0, 0, 1; t4, -t3, t14, t36; t2, -t1, t12, t35; t8, -t48, -t22, t38; 0, 0, 0, 1; t4, t14, t3, t4 * pkin(5) + t3 * qJ(6) + t36; t2, t12, t1, t2 * pkin(5) + t1 * qJ(6) + t35; t8, -t22, t48 (pkin(5) * t21 + qJ(6) * t19) * t20 + t38; 0, 0, 0, 1;];
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
