% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRPRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:06:34
% EndTime: 2018-11-23 15:06:35
% DurationCPUTime: 0.19s
% Computational Cost: add. (520->80), mult. (558->86), div. (0->0), fcn. (622->16), ass. (0->53)
t33 = sin(pkin(11));
t34 = sin(pkin(6));
t22 = t33 * t34;
t35 = cos(pkin(11));
t67 = t35 * t34;
t66 = pkin(6) - qJ(2);
t65 = pkin(6) + qJ(2);
t64 = pkin(7) * t67;
t63 = t33 * pkin(1) + 0;
t62 = qJ(1) + 0;
t37 = sin(qJ(5));
t61 = pkin(5) * t37 + pkin(8);
t60 = t35 * pkin(1) + pkin(7) * t22 + 0;
t36 = cos(pkin(6));
t59 = t36 * pkin(7) + t62;
t58 = cos(t65);
t57 = sin(t66);
t56 = cos(t66) / 0.2e1;
t55 = sin(t65) / 0.2e1;
t54 = t56 + t58 / 0.2e1;
t39 = sin(qJ(2));
t8 = t33 * t39 - t35 * t54;
t42 = cos(qJ(2));
t48 = t55 - t57 / 0.2e1;
t9 = t33 * t42 + t35 * t48;
t53 = t9 * pkin(2) + t8 * qJ(3) + t63;
t10 = t33 * t54 + t35 * t39;
t11 = -t33 * t48 + t35 * t42;
t52 = t11 * pkin(2) + t10 * qJ(3) + t60;
t18 = t55 + t57 / 0.2e1;
t19 = t56 - t58 / 0.2e1;
t51 = t19 * pkin(2) - t18 * qJ(3) + t59;
t50 = pkin(3) * t22 + t52;
t49 = t36 * pkin(3) + t51;
t47 = t11 * pkin(8) + t50;
t46 = (-pkin(3) - pkin(7)) * t67 + t53;
t45 = t19 * pkin(8) + t49;
t44 = t9 * pkin(8) + t46;
t43 = -pkin(10) - pkin(9);
t41 = cos(qJ(4));
t40 = cos(qJ(5));
t38 = sin(qJ(4));
t32 = qJ(5) + qJ(6);
t27 = cos(t32);
t26 = sin(t32);
t25 = t40 * pkin(5) + pkin(4);
t13 = -t18 * t38 + t36 * t41;
t12 = t18 * t41 + t36 * t38;
t4 = t8 * t38 - t41 * t67;
t3 = t38 * t67 + t8 * t41;
t2 = t10 * t38 + t41 * t22;
t1 = -t10 * t41 + t38 * t22;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t33, 0, 0; t33, t35, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t11, -t10, t22, t60; t9, -t8, -t67, t63 - t64; t19, t18, t36, t59; 0, 0, 0, 1; t22, -t11, t10, t52; -t67, -t9, t8, t53 - t64; t36, -t19, -t18, t51; 0, 0, 0, 1; t2, -t1, t11, t47; t4, t3, t9, t44; t13, -t12, t19, t45; 0, 0, 0, 1; t11 * t37 + t2 * t40, t11 * t40 - t2 * t37, t1, t2 * pkin(4) + t1 * pkin(9) + t47; t9 * t37 + t4 * t40, -t4 * t37 + t9 * t40, -t3, t4 * pkin(4) - t3 * pkin(9) + t44; t13 * t40 + t19 * t37, -t13 * t37 + t19 * t40, t12, t13 * pkin(4) + t12 * pkin(9) + t45; 0, 0, 0, 1; t11 * t26 + t2 * t27, t11 * t27 - t2 * t26, t1, -t1 * t43 + t61 * t11 + t2 * t25 + t50; t9 * t26 + t4 * t27, -t4 * t26 + t9 * t27, -t3, t4 * t25 + t3 * t43 + t61 * t9 + t46; t13 * t27 + t19 * t26, -t13 * t26 + t19 * t27, t12, -t12 * t43 + t13 * t25 + t61 * t19 + t49; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
