% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2018-11-23 15:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPRRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:45:11
% EndTime: 2018-11-23 15:45:11
% DurationCPUTime: 0.13s
% Computational Cost: add. (147->50), mult. (232->43), div. (0->0), fcn. (340->8), ass. (0->34)
t50 = sin(qJ(1));
t32 = cos(qJ(1));
t42 = sin(pkin(9));
t43 = cos(pkin(9));
t15 = -t32 * t43 - t50 * t42;
t29 = sin(qJ(4));
t11 = t15 * t29;
t31 = cos(qJ(4));
t49 = t15 * t31;
t16 = t32 * t42 - t50 * t43;
t10 = t16 * t29;
t48 = t16 * t31;
t28 = sin(qJ(5));
t47 = t28 * t31;
t46 = t29 * t28;
t30 = cos(qJ(5));
t45 = t29 * t30;
t44 = t30 * t31;
t27 = pkin(6) + 0;
t20 = -qJ(3) + t27;
t41 = t32 * pkin(1) + t50 * qJ(2) + 0;
t40 = t31 * pkin(8) + t20;
t39 = t32 * pkin(2) + t41;
t38 = t50 * pkin(1) - t32 * qJ(2) + 0;
t37 = t50 * pkin(2) + t38;
t36 = -t15 * pkin(3) + t16 * pkin(7) + t39;
t35 = -pkin(4) * t49 - pkin(8) * t11 + t36;
t34 = -t16 * pkin(3) - t15 * pkin(7) + t37;
t33 = -pkin(4) * t48 - pkin(8) * t10 + t34;
t4 = -t15 * t44 + t16 * t28;
t3 = -t15 * t47 - t16 * t30;
t2 = -t15 * t28 - t16 * t44;
t1 = t15 * t30 - t16 * t47;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t32, -t50, 0, 0; t50, t32, 0, 0; 0, 0, 1, t27; 0, 0, 0, 1; t32, 0, t50, t41; t50, 0, -t32, t38; 0, 1, 0, t27; 0, 0, 0, 1; -t15, -t16, 0, t39; -t16, t15, 0, t37; 0, 0, -1, t20; 0, 0, 0, 1; -t49, t11, t16, t36; -t48, t10, -t15, t34; -t29, -t31, 0, t20; 0, 0, 0, 1; t4, -t3, -t11, t35; t2, -t1, -t10, t33; -t45, t46, t31, -t29 * pkin(4) + t40; 0, 0, 0, 1; t4, -t11, t3, t4 * pkin(5) + t3 * qJ(6) + t35; t2, -t10, t1, t2 * pkin(5) + t1 * qJ(6) + t33; -t45, t31, -t46 (-pkin(5) * t30 - qJ(6) * t28 - pkin(4)) * t29 + t40; 0, 0, 0, 1;];
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
