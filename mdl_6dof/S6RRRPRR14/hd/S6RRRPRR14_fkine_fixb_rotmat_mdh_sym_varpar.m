% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 18:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPRR14_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:01:56
% EndTime: 2018-11-23 18:01:56
% DurationCPUTime: 0.20s
% Computational Cost: add. (584->81), mult. (646->87), div. (0->0), fcn. (727->16), ass. (0->57)
t73 = pkin(4) + pkin(9);
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t64 = pkin(6) - qJ(2);
t53 = cos(t64) / 0.2e1;
t63 = pkin(6) + qJ(2);
t57 = cos(t63);
t44 = t53 + t57 / 0.2e1;
t12 = t38 * t37 - t41 * t44;
t72 = t12 * pkin(9);
t14 = t41 * t37 + t38 * t44;
t71 = t14 * pkin(9);
t52 = sin(t63) / 0.2e1;
t56 = sin(t64);
t18 = t52 + t56 / 0.2e1;
t70 = t18 * pkin(9);
t39 = cos(qJ(5));
t69 = t39 * pkin(5) + t73;
t68 = cos(qJ(3));
t34 = sin(pkin(6));
t67 = t38 * t34;
t66 = t41 * t34;
t65 = cos(pkin(6));
t62 = pkin(7) + 0;
t61 = t34 * t68;
t60 = t65 * pkin(8) + t62;
t35 = sin(qJ(5));
t59 = pkin(5) * t35 + qJ(4);
t58 = t41 * pkin(1) + pkin(8) * t67 + 0;
t20 = t53 - t57 / 0.2e1;
t55 = t20 * pkin(2) + t60;
t19 = t52 - t56 / 0.2e1;
t40 = cos(qJ(2));
t15 = -t38 * t19 + t41 * t40;
t54 = t15 * pkin(2) + t58;
t36 = sin(qJ(3));
t11 = t20 * t68 + t65 * t36;
t51 = t11 * pkin(3) + t55;
t6 = t15 * t68 + t36 * t67;
t50 = t6 * pkin(3) + t54;
t49 = t38 * pkin(1) - pkin(8) * t66 + 0;
t13 = t41 * t19 + t38 * t40;
t48 = t13 * pkin(2) + t49;
t4 = t13 * t68 - t36 * t66;
t47 = t4 * pkin(3) + t48;
t5 = t15 * t36 - t38 * t61;
t46 = t5 * qJ(4) + t50;
t10 = t20 * t36 - t65 * t68;
t45 = t10 * qJ(4) + t51;
t3 = t13 * t36 + t41 * t61;
t43 = t3 * qJ(4) + t47;
t42 = -pkin(11) - pkin(10);
t33 = qJ(5) + qJ(6);
t29 = cos(t33);
t28 = sin(t33);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t41, -t38, 0, 0; t38, t41, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t15, -t14, t67, t58; t13, -t12, -t66, t49; t20, t18, t65, t60; 0, 0, 0, 1; t6, -t5, t14, t54 + t71; t4, -t3, t12, t48 + t72; t11, -t10, -t18, t55 - t70; 0, 0, 0, 1; t14, -t6, t5, t46 + t71; t12, -t4, t3, t43 + t72; -t18, -t11, t10, t45 - t70; 0, 0, 0, 1; t14 * t39 + t5 * t35, -t14 * t35 + t5 * t39, t6, t6 * pkin(10) + t73 * t14 + t46; t12 * t39 + t3 * t35, -t12 * t35 + t3 * t39, t4, t4 * pkin(10) + t73 * t12 + t43; t10 * t35 - t18 * t39, t10 * t39 + t18 * t35, t11, t11 * pkin(10) - t73 * t18 + t45; 0, 0, 0, 1; t14 * t29 + t5 * t28, -t14 * t28 + t5 * t29, t6, t69 * t14 - t6 * t42 + t59 * t5 + t50; t12 * t29 + t3 * t28, -t12 * t28 + t3 * t29, t4, t69 * t12 + t59 * t3 - t4 * t42 + t47; t10 * t28 - t18 * t29, t10 * t29 + t18 * t28, t11, t59 * t10 - t11 * t42 - t69 * t18 + t51; 0, 0, 0, 1;];
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
