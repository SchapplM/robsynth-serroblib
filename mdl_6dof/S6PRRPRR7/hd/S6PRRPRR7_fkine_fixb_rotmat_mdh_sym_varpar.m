% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRRPRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:18:37
% EndTime: 2018-11-23 15:18:37
% DurationCPUTime: 0.22s
% Computational Cost: add. (584->81), mult. (646->87), div. (0->0), fcn. (727->16), ass. (0->57)
t73 = pkin(4) + pkin(8);
t34 = sin(pkin(11));
t36 = cos(pkin(11));
t39 = sin(qJ(2));
t64 = pkin(6) - qJ(2);
t53 = cos(t64) / 0.2e1;
t63 = pkin(6) + qJ(2);
t57 = cos(t63);
t44 = t53 + t57 / 0.2e1;
t10 = t34 * t39 - t36 * t44;
t72 = t10 * pkin(8);
t12 = t34 * t44 + t36 * t39;
t71 = t12 * pkin(8);
t52 = sin(t63) / 0.2e1;
t56 = sin(t64);
t18 = t52 + t56 / 0.2e1;
t70 = t18 * pkin(8);
t40 = cos(qJ(5));
t69 = t40 * pkin(5) + t73;
t68 = cos(qJ(3));
t35 = sin(pkin(6));
t67 = t34 * t35;
t66 = t36 * t35;
t65 = cos(pkin(6));
t62 = qJ(1) + 0;
t61 = t35 * t68;
t37 = sin(qJ(5));
t60 = pkin(5) * t37 + qJ(4);
t59 = t36 * pkin(1) + pkin(7) * t67 + 0;
t58 = t65 * pkin(7) + t62;
t19 = t52 - t56 / 0.2e1;
t41 = cos(qJ(2));
t13 = -t34 * t19 + t36 * t41;
t55 = t13 * pkin(2) + t59;
t20 = t53 - t57 / 0.2e1;
t54 = t20 * pkin(2) + t58;
t38 = sin(qJ(3));
t6 = t13 * t68 + t38 * t67;
t51 = t6 * pkin(3) + t55;
t15 = t20 * t68 + t65 * t38;
t50 = t15 * pkin(3) + t54;
t49 = t34 * pkin(1) - pkin(7) * t66 + 0;
t11 = t36 * t19 + t34 * t41;
t48 = t11 * pkin(2) + t49;
t4 = t11 * t68 - t38 * t66;
t47 = t4 * pkin(3) + t48;
t5 = t13 * t38 - t34 * t61;
t46 = t5 * qJ(4) + t51;
t14 = t20 * t38 - t65 * t68;
t45 = t14 * qJ(4) + t50;
t3 = t11 * t38 + t36 * t61;
t43 = t3 * qJ(4) + t47;
t42 = -pkin(10) - pkin(9);
t33 = qJ(5) + qJ(6);
t29 = cos(t33);
t28 = sin(t33);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t36, -t34, 0, 0; t34, t36, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t13, -t12, t67, t59; t11, -t10, -t66, t49; t20, t18, t65, t58; 0, 0, 0, 1; t6, -t5, t12, t55 + t71; t4, -t3, t10, t48 + t72; t15, -t14, -t18, t54 - t70; 0, 0, 0, 1; t12, -t6, t5, t46 + t71; t10, -t4, t3, t43 + t72; -t18, -t15, t14, t45 - t70; 0, 0, 0, 1; t12 * t40 + t5 * t37, -t12 * t37 + t5 * t40, t6, t6 * pkin(9) + t73 * t12 + t46; t10 * t40 + t3 * t37, -t10 * t37 + t3 * t40, t4, t4 * pkin(9) + t73 * t10 + t43; t14 * t37 - t18 * t40, t14 * t40 + t18 * t37, t15, t15 * pkin(9) - t73 * t18 + t45; 0, 0, 0, 1; t12 * t29 + t5 * t28, -t12 * t28 + t5 * t29, t6, t69 * t12 - t6 * t42 + t60 * t5 + t51; t10 * t29 + t3 * t28, -t10 * t28 + t3 * t29, t4, t69 * t10 + t60 * t3 - t4 * t42 + t47; t14 * t28 - t18 * t29, t14 * t29 + t18 * t28, t15, t60 * t14 - t15 * t42 - t69 * t18 + t50; 0, 0, 0, 1;];
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
