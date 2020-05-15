% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRRR13_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:29:48
% EndTime: 2018-11-23 17:29:49
% DurationCPUTime: 0.18s
% Computational Cost: add. (520->80), mult. (558->86), div. (0->0), fcn. (622->16), ass. (0->53)
t33 = sin(pkin(6));
t38 = sin(qJ(1));
t24 = t38 * t33;
t42 = cos(qJ(1));
t67 = t42 * t33;
t66 = pkin(6) - qJ(2);
t65 = pkin(6) + qJ(2);
t64 = pkin(7) + 0;
t63 = pkin(8) * t67;
t62 = t38 * pkin(1) + 0;
t35 = sin(qJ(5));
t61 = pkin(5) * t35 + pkin(9);
t34 = cos(pkin(6));
t60 = t34 * pkin(8) + t64;
t59 = t42 * pkin(1) + pkin(8) * t24 + 0;
t58 = cos(t65);
t57 = sin(t66);
t56 = cos(t66) / 0.2e1;
t55 = sin(t65) / 0.2e1;
t54 = t56 + t58 / 0.2e1;
t37 = sin(qJ(2));
t10 = t38 * t37 - t42 * t54;
t41 = cos(qJ(2));
t48 = t55 - t57 / 0.2e1;
t11 = t38 * t41 + t42 * t48;
t53 = t11 * pkin(2) + t10 * qJ(3) + t62;
t18 = t55 + t57 / 0.2e1;
t19 = t56 - t58 / 0.2e1;
t52 = t19 * pkin(2) - t18 * qJ(3) + t60;
t12 = t42 * t37 + t38 * t54;
t13 = -t38 * t48 + t42 * t41;
t51 = t13 * pkin(2) + t12 * qJ(3) + t59;
t50 = t34 * pkin(3) + t52;
t49 = pkin(3) * t24 + t51;
t47 = t19 * pkin(9) + t50;
t46 = t13 * pkin(9) + t49;
t45 = (-pkin(3) - pkin(8)) * t67 + t53;
t44 = t11 * pkin(9) + t45;
t43 = -pkin(11) - pkin(10);
t40 = cos(qJ(4));
t39 = cos(qJ(5));
t36 = sin(qJ(4));
t32 = qJ(5) + qJ(6);
t27 = cos(t32);
t26 = sin(t32);
t25 = t39 * pkin(5) + pkin(4);
t9 = -t18 * t36 + t34 * t40;
t8 = t18 * t40 + t34 * t36;
t4 = t10 * t36 - t40 * t67;
t3 = t10 * t40 + t36 * t67;
t2 = t12 * t36 + t40 * t24;
t1 = -t12 * t40 + t36 * t24;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t42, -t38, 0, 0; t38, t42, 0, 0; 0, 0, 1, t64; 0, 0, 0, 1; t13, -t12, t24, t59; t11, -t10, -t67, t62 - t63; t19, t18, t34, t60; 0, 0, 0, 1; t24, -t13, t12, t51; -t67, -t11, t10, t53 - t63; t34, -t19, -t18, t52; 0, 0, 0, 1; t2, -t1, t13, t46; t4, t3, t11, t44; t9, -t8, t19, t47; 0, 0, 0, 1; t13 * t35 + t2 * t39, t13 * t39 - t2 * t35, t1, t2 * pkin(4) + t1 * pkin(10) + t46; t11 * t35 + t4 * t39, t11 * t39 - t4 * t35, -t3, t4 * pkin(4) - t3 * pkin(10) + t44; t19 * t35 + t9 * t39, t19 * t39 - t9 * t35, t8, t9 * pkin(4) + t8 * pkin(10) + t47; 0, 0, 0, 1; t13 * t26 + t2 * t27, t13 * t27 - t2 * t26, t1, -t1 * t43 + t61 * t13 + t2 * t25 + t49; t11 * t26 + t4 * t27, t11 * t27 - t4 * t26, -t3, t61 * t11 + t4 * t25 + t3 * t43 + t45; t19 * t26 + t9 * t27, t19 * t27 - t9 * t26, t8, t61 * t19 + t9 * t25 - t8 * t43 + t50; 0, 0, 0, 1;];
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
