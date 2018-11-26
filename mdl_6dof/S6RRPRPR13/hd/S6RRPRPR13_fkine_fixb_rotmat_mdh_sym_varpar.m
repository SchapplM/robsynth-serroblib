% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRPR13_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:09:06
% EndTime: 2018-11-23 17:09:06
% DurationCPUTime: 0.18s
% Computational Cost: add. (520->80), mult. (558->86), div. (0->0), fcn. (622->16), ass. (0->53)
t34 = sin(pkin(6));
t40 = sin(qJ(1));
t24 = t40 * t34;
t43 = cos(qJ(1));
t67 = t43 * t34;
t66 = pkin(6) - qJ(2);
t65 = pkin(6) + qJ(2);
t64 = pkin(7) + 0;
t63 = pkin(8) * t67;
t62 = t40 * pkin(1) + 0;
t33 = sin(pkin(11));
t61 = pkin(5) * t33 + pkin(9);
t36 = cos(pkin(6));
t60 = t36 * pkin(8) + t64;
t59 = t43 * pkin(1) + pkin(8) * t24 + 0;
t58 = cos(t65);
t57 = sin(t66);
t56 = cos(t66) / 0.2e1;
t55 = sin(t65) / 0.2e1;
t54 = t56 + t58 / 0.2e1;
t39 = sin(qJ(2));
t10 = t40 * t39 - t43 * t54;
t42 = cos(qJ(2));
t48 = t55 - t57 / 0.2e1;
t11 = t40 * t42 + t43 * t48;
t53 = t11 * pkin(2) + t10 * qJ(3) + t62;
t18 = t55 + t57 / 0.2e1;
t19 = t56 - t58 / 0.2e1;
t52 = t19 * pkin(2) - t18 * qJ(3) + t60;
t12 = t43 * t39 + t40 * t54;
t13 = -t40 * t48 + t43 * t42;
t51 = t13 * pkin(2) + t12 * qJ(3) + t59;
t50 = t36 * pkin(3) + t52;
t49 = pkin(3) * t24 + t51;
t47 = t19 * pkin(9) + t50;
t46 = t13 * pkin(9) + t49;
t45 = (-pkin(3) - pkin(8)) * t67 + t53;
t44 = t11 * pkin(9) + t45;
t41 = cos(qJ(4));
t38 = sin(qJ(4));
t37 = -pkin(10) - qJ(5);
t35 = cos(pkin(11));
t32 = pkin(11) + qJ(6);
t27 = cos(t32);
t26 = sin(t32);
t25 = t35 * pkin(5) + pkin(4);
t9 = -t18 * t38 + t36 * t41;
t8 = t18 * t41 + t36 * t38;
t4 = t10 * t38 - t41 * t67;
t3 = t10 * t41 + t38 * t67;
t2 = t12 * t38 + t41 * t24;
t1 = -t12 * t41 + t38 * t24;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t43, -t40, 0, 0; t40, t43, 0, 0; 0, 0, 1, t64; 0, 0, 0, 1; t13, -t12, t24, t59; t11, -t10, -t67, t62 - t63; t19, t18, t36, t60; 0, 0, 0, 1; t24, -t13, t12, t51; -t67, -t11, t10, t53 - t63; t36, -t19, -t18, t52; 0, 0, 0, 1; t2, -t1, t13, t46; t4, t3, t11, t44; t9, -t8, t19, t47; 0, 0, 0, 1; t13 * t33 + t2 * t35, t13 * t35 - t2 * t33, t1, t2 * pkin(4) + t1 * qJ(5) + t46; t11 * t33 + t4 * t35, t11 * t35 - t4 * t33, -t3, t4 * pkin(4) - t3 * qJ(5) + t44; t19 * t33 + t9 * t35, t19 * t35 - t9 * t33, t8, t9 * pkin(4) + t8 * qJ(5) + t47; 0, 0, 0, 1; t13 * t26 + t2 * t27, t13 * t27 - t2 * t26, t1, -t1 * t37 + t61 * t13 + t2 * t25 + t49; t11 * t26 + t4 * t27, t11 * t27 - t4 * t26, -t3, t61 * t11 + t4 * t25 + t3 * t37 + t45; t19 * t26 + t9 * t27, t19 * t27 - t9 * t26, t8, t61 * t19 + t9 * t25 - t8 * t37 + t50; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
