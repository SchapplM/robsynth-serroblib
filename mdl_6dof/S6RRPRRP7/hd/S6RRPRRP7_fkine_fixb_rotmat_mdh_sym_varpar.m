% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRRP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:15:01
% EndTime: 2018-11-23 17:15:02
% DurationCPUTime: 0.13s
% Computational Cost: add. (145->49), mult. (277->53), div. (0->0), fcn. (377->8), ass. (0->39)
t35 = sin(qJ(2));
t38 = cos(qJ(4));
t34 = sin(qJ(4));
t39 = cos(qJ(2));
t54 = t39 * t34;
t14 = t35 * t38 - t54;
t33 = sin(qJ(5));
t56 = t14 * t33;
t36 = sin(qJ(1));
t55 = t36 * t35;
t23 = t36 * t39;
t40 = cos(qJ(1));
t53 = t40 * t35;
t25 = t40 * t39;
t52 = qJ(3) * t35;
t32 = pkin(6) + 0;
t51 = t40 * pkin(1) + t36 * pkin(7) + 0;
t13 = t35 * t34 + t39 * t38;
t50 = t36 * pkin(1) - t40 * pkin(7) + 0;
t49 = pkin(2) * t25 + t40 * t52 + t51;
t48 = t35 * pkin(2) - t39 * qJ(3) + t32;
t47 = pkin(2) * t23 + t36 * t52 + t50;
t46 = t35 * pkin(3) + t48;
t45 = pkin(3) * t23 + t40 * pkin(8) + t47;
t44 = pkin(3) * t25 - t36 * pkin(8) + t49;
t43 = t14 * pkin(4) + t13 * pkin(9) + t46;
t7 = t36 * t54 - t38 * t55;
t8 = t13 * t36;
t42 = t8 * pkin(4) + t7 * pkin(9) + t45;
t10 = t13 * t40;
t9 = t34 * t25 - t38 * t53;
t41 = t10 * pkin(4) + t9 * pkin(9) + t44;
t37 = cos(qJ(5));
t11 = t14 * t37;
t4 = t10 * t37 - t36 * t33;
t3 = t10 * t33 + t36 * t37;
t2 = t40 * t33 + t8 * t37;
t1 = t8 * t33 - t40 * t37;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t40, -t36, 0, 0; t36, t40, 0, 0; 0, 0, 1, t32; 0, 0, 0, 1; t25, -t53, t36, t51; t23, -t55, -t40, t50; t35, t39, 0, t32; 0, 0, 0, 1; t25, t36, t53, t49; t23, -t40, t55, t47; t35, 0, -t39, t48; 0, 0, 0, 1; t10, -t9, -t36, t44; t8, -t7, t40, t45; t14, -t13, 0, t46; 0, 0, 0, 1; t4, -t3, t9, t41; t2, -t1, t7, t42; t11, -t56, t13, t43; 0, 0, 0, 1; t4, t9, t3, t4 * pkin(5) + t3 * qJ(6) + t41; t2, t7, t1, t2 * pkin(5) + t1 * qJ(6) + t42; t11, t13, t56 (pkin(5) * t37 + qJ(6) * t33) * t14 + t43; 0, 0, 0, 1;];
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
