% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2018-11-23 16:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPPRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:47:18
% EndTime: 2018-11-23 16:47:19
% DurationCPUTime: 0.14s
% Computational Cost: add. (147->62), mult. (177->58), div. (0->0), fcn. (240->8), ass. (0->39)
t28 = sin(qJ(1));
t29 = cos(qJ(2));
t14 = t28 * t29;
t27 = sin(qJ(2));
t45 = qJ(3) * t27;
t50 = pkin(2) * t14 + t28 * t45;
t49 = t28 * t27;
t22 = pkin(9) + qJ(5);
t16 = sin(t22);
t48 = t29 * t16;
t17 = cos(t22);
t47 = t29 * t17;
t30 = cos(qJ(1));
t46 = t30 * t27;
t15 = t30 * t29;
t44 = qJ(4) * t29;
t23 = pkin(6) + 0;
t43 = t28 * pkin(1) + 0;
t24 = sin(pkin(9));
t42 = t24 * t49;
t41 = t24 * t46;
t40 = t27 * pkin(2) + t23;
t39 = -pkin(4) * t24 - qJ(3);
t38 = t30 * pkin(1) + t28 * pkin(7) + 0;
t37 = t43 + t50;
t36 = -t30 * pkin(7) + t43;
t35 = pkin(2) * t15 + t30 * t45 + t38;
t26 = -pkin(8) - qJ(4);
t34 = -t27 * t26 + t40;
t33 = -t29 * qJ(3) + t40;
t25 = cos(pkin(9));
t13 = t25 * pkin(4) + pkin(3);
t32 = pkin(4) * t41 + t28 * t13 - t26 * t15 + t35;
t31 = -t26 * t14 + pkin(4) * t42 + (-pkin(7) - t13) * t30 + t37;
t4 = t16 * t49 - t30 * t17;
t3 = t30 * t16 + t17 * t49;
t2 = t16 * t46 + t28 * t17;
t1 = t28 * t16 - t17 * t46;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t28, 0, 0; t28, t30, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t15, -t46, t28, t38; t14, -t49, -t30, t36; t27, t29, 0, t23; 0, 0, 0, 1; t28, -t15, t46, t35; -t30, -t14, t49, t36 + t50; 0, -t27, -t29, t33; 0, 0, 0, 1; t28 * t25 + t41, -t28 * t24 + t25 * t46, t15, t28 * pkin(3) + t30 * t44 + t35; -t30 * t25 + t42, t30 * t24 + t25 * t49, t14, t28 * t44 + (-pkin(3) - pkin(7)) * t30 + t37; -t29 * t24, -t29 * t25, t27, t27 * qJ(4) + t33; 0, 0, 0, 1; t2, -t1, t15, t32; t4, t3, t14, t31; -t48, -t47, t27, t39 * t29 + t34; 0, 0, 0, 1; t2, t15, t1, t2 * pkin(5) + t1 * qJ(6) + t32; t4, t14, -t3, t4 * pkin(5) - t3 * qJ(6) + t31; -t48, t27, t47 (-pkin(5) * t16 + qJ(6) * t17 + t39) * t29 + t34; 0, 0, 0, 1;];
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
