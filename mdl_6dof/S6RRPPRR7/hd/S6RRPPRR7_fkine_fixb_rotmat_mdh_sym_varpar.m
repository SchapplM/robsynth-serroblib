% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2018-11-23 16:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPPRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:53:05
% EndTime: 2018-11-23 16:53:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (464->67), mult. (485->70), div. (0->0), fcn. (530->14), ass. (0->47)
t34 = sin(pkin(6));
t39 = sin(qJ(1));
t30 = t39 * t34;
t43 = cos(qJ(1));
t66 = t43 * t34;
t65 = qJ(4) * t34;
t64 = pkin(6) - qJ(2);
t63 = pkin(6) + qJ(2);
t62 = pkin(7) + 0;
t35 = cos(pkin(6));
t61 = t35 * pkin(8) + t62;
t60 = t43 * pkin(1) + pkin(8) * t30 + 0;
t59 = cos(t63);
t58 = sin(t64);
t57 = cos(t64) / 0.2e1;
t56 = sin(t63) / 0.2e1;
t55 = t56 - t58 / 0.2e1;
t54 = t39 * pkin(1) - pkin(8) * t66 + 0;
t22 = t56 + t58 / 0.2e1;
t23 = t57 - t59 / 0.2e1;
t53 = t23 * pkin(2) - t22 * qJ(3) + t61;
t38 = sin(qJ(2));
t50 = t57 + t59 / 0.2e1;
t15 = t43 * t38 + t39 * t50;
t42 = cos(qJ(2));
t16 = -t39 * t55 + t43 * t42;
t52 = t16 * pkin(2) + t15 * qJ(3) + t60;
t13 = t39 * t38 - t43 * t50;
t14 = t39 * t42 + t43 * t55;
t51 = t14 * pkin(2) + t13 * qJ(3) + t54;
t49 = t23 * pkin(3) - t35 * qJ(4) + t53;
t48 = t14 * pkin(3) + t43 * t65 + t51;
t47 = t16 * pkin(3) - t39 * t65 + t52;
t46 = -t22 * pkin(4) + t23 * pkin(9) + t49;
t45 = t13 * pkin(4) + t14 * pkin(9) + t48;
t44 = t15 * pkin(4) + t16 * pkin(9) + t47;
t41 = cos(qJ(5));
t40 = cos(qJ(6));
t37 = sin(qJ(5));
t36 = sin(qJ(6));
t12 = -t22 * t41 - t35 * t37;
t11 = -t22 * t37 + t35 * t41;
t4 = t15 * t41 - t37 * t30;
t3 = t15 * t37 + t41 * t30;
t2 = t13 * t41 + t37 * t66;
t1 = t13 * t37 - t41 * t66;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t43, -t39, 0, 0; t39, t43, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t16, -t15, t30, t60; t14, -t13, -t66, t54; t23, t22, t35, t61; 0, 0, 0, 1; t16, t30, t15, t52; t14, -t66, t13, t51; t23, t35, -t22, t53; 0, 0, 0, 1; t15, -t16, -t30, t47; t13, -t14, t66, t48; -t22, -t23, -t35, t49; 0, 0, 0, 1; t4, -t3, t16, t44; t2, -t1, t14, t45; t12, -t11, t23, t46; 0, 0, 0, 1; t16 * t36 + t4 * t40, t16 * t40 - t4 * t36, t3, t4 * pkin(5) + t3 * pkin(10) + t44; t14 * t36 + t2 * t40, t14 * t40 - t2 * t36, t1, t2 * pkin(5) + t1 * pkin(10) + t45; t12 * t40 + t23 * t36, -t12 * t36 + t23 * t40, t11, t12 * pkin(5) + t11 * pkin(10) + t46; 0, 0, 0, 1;];
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
