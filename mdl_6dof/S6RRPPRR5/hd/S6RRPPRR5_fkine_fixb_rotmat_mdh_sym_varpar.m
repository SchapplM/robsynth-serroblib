% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR5
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
% Datum: 2018-11-23 16:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPPRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:51:31
% EndTime: 2018-11-23 16:51:31
% DurationCPUTime: 0.18s
% Computational Cost: add. (464->70), mult. (485->70), div. (0->0), fcn. (530->14), ass. (0->51)
t61 = pkin(6) + qJ(2);
t52 = sin(t61) / 0.2e1;
t62 = pkin(6) - qJ(2);
t56 = sin(t62);
t22 = t52 - t56 / 0.2e1;
t38 = sin(qJ(1));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t14 = t42 * t22 + t38 * t41;
t33 = sin(pkin(6));
t65 = qJ(4) * t33;
t68 = t14 * pkin(3) + t42 * t65;
t29 = t38 * t33;
t67 = t42 * t33;
t66 = pkin(9) - qJ(3);
t37 = sin(qJ(2));
t53 = cos(t62) / 0.2e1;
t57 = cos(t61);
t46 = t53 + t57 / 0.2e1;
t15 = t42 * t37 + t38 * t46;
t64 = t15 * qJ(3);
t21 = t52 + t56 / 0.2e1;
t63 = t21 * qJ(3);
t60 = pkin(7) + 0;
t34 = cos(pkin(6));
t59 = t34 * pkin(8) + t60;
t58 = t42 * pkin(1) + pkin(8) * t29 + 0;
t23 = t53 - t57 / 0.2e1;
t55 = t23 * pkin(2) + t59;
t16 = -t38 * t22 + t42 * t41;
t54 = t16 * pkin(2) + t58;
t51 = t38 * pkin(1) - pkin(8) * t67 + 0;
t50 = t14 * pkin(2) + t51;
t49 = t23 * pkin(3) - t34 * qJ(4) + t55;
t48 = t16 * pkin(3) - t38 * t65 + t54;
t13 = t38 * t37 - t42 * t46;
t47 = t13 * qJ(3) + t50;
t45 = t23 * pkin(4) + t66 * t21 + t49;
t44 = t14 * pkin(4) - t66 * t13 + t50 + t68;
t43 = t16 * pkin(4) - t66 * t15 + t48;
t40 = cos(qJ(5));
t39 = cos(qJ(6));
t36 = sin(qJ(5));
t35 = sin(qJ(6));
t12 = t23 * t40 - t34 * t36;
t11 = t23 * t36 + t34 * t40;
t4 = t16 * t40 - t36 * t29;
t3 = t16 * t36 + t40 * t29;
t2 = t14 * t40 + t36 * t67;
t1 = t14 * t36 - t40 * t67;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t42, -t38, 0, 0; t38, t42, 0, 0; 0, 0, 1, t60; 0, 0, 0, 1; t16, -t15, t29, t58; t14, -t13, -t67, t51; t23, t21, t34, t59; 0, 0, 0, 1; t16, t29, t15, t54 + t64; t14, -t67, t13, t47; t23, t34, -t21, t55 - t63; 0, 0, 0, 1; t16, t15, -t29, t48 + t64; t14, t13, t67, t47 + t68; t23, -t21, -t34, t49 - t63; 0, 0, 0, 1; t4, -t3, -t15, t43; t2, -t1, -t13, t44; t12, -t11, t21, t45; 0, 0, 0, 1; -t15 * t35 + t4 * t39, -t15 * t39 - t4 * t35, t3, t4 * pkin(5) + t3 * pkin(10) + t43; -t13 * t35 + t2 * t39, -t13 * t39 - t2 * t35, t1, t2 * pkin(5) + t1 * pkin(10) + t44; t12 * t39 + t21 * t35, -t12 * t35 + t21 * t39, t11, t12 * pkin(5) + t11 * pkin(10) + t45; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
