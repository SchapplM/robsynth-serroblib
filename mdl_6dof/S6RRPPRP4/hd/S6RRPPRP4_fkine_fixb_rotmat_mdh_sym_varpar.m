% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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

function [T_c_mdh, Tc_stack] = S6RRPPRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:46:35
% EndTime: 2018-11-23 16:46:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (163->55), mult. (316->56), div. (0->0), fcn. (430->8), ass. (0->38)
t35 = sin(pkin(9));
t38 = sin(qJ(2));
t57 = t38 * t35;
t36 = cos(pkin(9));
t25 = t38 * t36;
t39 = sin(qJ(1));
t27 = t39 * t38;
t41 = cos(qJ(2));
t56 = t39 * t41;
t42 = cos(qJ(1));
t28 = t42 * t38;
t55 = t42 * t41;
t54 = qJ(3) * t38;
t34 = pkin(6) + 0;
t53 = t42 * pkin(1) + t39 * pkin(7) + 0;
t52 = t39 * pkin(1) - t42 * pkin(7) + 0;
t51 = pkin(2) * t55 + t42 * t54 + t53;
t50 = t38 * pkin(2) - t41 * qJ(3) + t34;
t49 = pkin(2) * t56 + t39 * t54 + t52;
t48 = pkin(3) * t25 + qJ(4) * t57 + t50;
t15 = t35 * t55 - t39 * t36;
t16 = t39 * t35 + t36 * t55;
t47 = t16 * pkin(3) + t15 * qJ(4) + t51;
t46 = pkin(4) * t25 + t41 * pkin(8) + t48;
t13 = t35 * t56 + t42 * t36;
t14 = -t42 * t35 + t36 * t56;
t45 = t14 * pkin(3) + t13 * qJ(4) + t49;
t44 = t16 * pkin(4) - pkin(8) * t28 + t47;
t43 = t14 * pkin(4) - pkin(8) * t27 + t45;
t40 = cos(qJ(5));
t37 = sin(qJ(5));
t8 = (t35 * t37 + t36 * t40) * t38;
t7 = t37 * t25 - t40 * t57;
t4 = t15 * t37 + t16 * t40;
t3 = -t15 * t40 + t16 * t37;
t2 = t13 * t37 + t14 * t40;
t1 = -t13 * t40 + t14 * t37;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t42, -t39, 0, 0; t39, t42, 0, 0; 0, 0, 1, t34; 0, 0, 0, 1; t55, -t28, t39, t53; t56, -t27, -t42, t52; t38, t41, 0, t34; 0, 0, 0, 1; t16, -t15, t28, t51; t14, -t13, t27, t49; t25, -t57, -t41, t50; 0, 0, 0, 1; t16, t28, t15, t47; t14, t27, t13, t45; t25, -t41, t57, t48; 0, 0, 0, 1; t4, -t3, -t28, t44; t2, -t1, -t27, t43; t8, -t7, t41, t46; 0, 0, 0, 1; t4, -t28, t3, t4 * pkin(5) + t3 * qJ(6) + t44; t2, -t27, t1, t2 * pkin(5) + t1 * qJ(6) + t43; t8, t41, t7, t8 * pkin(5) + t7 * qJ(6) + t46; 0, 0, 0, 1;];
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
