% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 14:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:58:14
% EndTime: 2018-11-23 14:58:14
% DurationCPUTime: 0.19s
% Computational Cost: add. (520->80), mult. (558->86), div. (0->0), fcn. (622->16), ass. (0->53)
t34 = sin(pkin(10));
t35 = sin(pkin(6));
t22 = t34 * t35;
t37 = cos(pkin(10));
t67 = t37 * t35;
t66 = pkin(6) - qJ(2);
t65 = pkin(6) + qJ(2);
t64 = pkin(7) * t67;
t63 = t34 * pkin(1) + 0;
t62 = qJ(1) + 0;
t33 = sin(pkin(11));
t61 = pkin(5) * t33 + pkin(8);
t60 = t37 * pkin(1) + pkin(7) * t22 + 0;
t38 = cos(pkin(6));
t59 = t38 * pkin(7) + t62;
t58 = cos(t65);
t57 = sin(t66);
t56 = cos(t66) / 0.2e1;
t55 = sin(t65) / 0.2e1;
t54 = t56 + t58 / 0.2e1;
t41 = sin(qJ(2));
t8 = t34 * t41 - t37 * t54;
t43 = cos(qJ(2));
t48 = t55 - t57 / 0.2e1;
t9 = t34 * t43 + t37 * t48;
t53 = t9 * pkin(2) + t8 * qJ(3) + t63;
t10 = t34 * t54 + t37 * t41;
t11 = -t34 * t48 + t37 * t43;
t52 = t11 * pkin(2) + t10 * qJ(3) + t60;
t18 = t55 + t57 / 0.2e1;
t19 = t56 - t58 / 0.2e1;
t51 = t19 * pkin(2) - t18 * qJ(3) + t59;
t50 = pkin(3) * t22 + t52;
t49 = t38 * pkin(3) + t51;
t47 = t11 * pkin(8) + t50;
t46 = (-pkin(3) - pkin(7)) * t67 + t53;
t45 = t19 * pkin(8) + t49;
t44 = t9 * pkin(8) + t46;
t42 = cos(qJ(4));
t40 = sin(qJ(4));
t39 = -pkin(9) - qJ(5);
t36 = cos(pkin(11));
t32 = pkin(11) + qJ(6);
t27 = cos(t32);
t26 = sin(t32);
t25 = t36 * pkin(5) + pkin(4);
t13 = -t18 * t40 + t38 * t42;
t12 = t18 * t42 + t38 * t40;
t4 = t8 * t40 - t42 * t67;
t3 = t40 * t67 + t8 * t42;
t2 = t10 * t40 + t42 * t22;
t1 = -t10 * t42 + t40 * t22;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t37, -t34, 0, 0; t34, t37, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t11, -t10, t22, t60; t9, -t8, -t67, t63 - t64; t19, t18, t38, t59; 0, 0, 0, 1; t22, -t11, t10, t52; -t67, -t9, t8, t53 - t64; t38, -t19, -t18, t51; 0, 0, 0, 1; t2, -t1, t11, t47; t4, t3, t9, t44; t13, -t12, t19, t45; 0, 0, 0, 1; t11 * t33 + t2 * t36, t11 * t36 - t2 * t33, t1, t2 * pkin(4) + t1 * qJ(5) + t47; t9 * t33 + t4 * t36, -t4 * t33 + t9 * t36, -t3, t4 * pkin(4) - t3 * qJ(5) + t44; t13 * t36 + t19 * t33, -t13 * t33 + t19 * t36, t12, t13 * pkin(4) + t12 * qJ(5) + t45; 0, 0, 0, 1; t11 * t26 + t2 * t27, t11 * t27 - t2 * t26, t1, -t1 * t39 + t61 * t11 + t2 * t25 + t50; t9 * t26 + t4 * t27, -t4 * t26 + t9 * t27, -t3, t4 * t25 + t3 * t39 + t61 * t9 + t46; t13 * t27 + t19 * t26, -t13 * t26 + t19 * t27, t12, -t12 * t39 + t13 * t25 + t61 * t19 + t49; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
