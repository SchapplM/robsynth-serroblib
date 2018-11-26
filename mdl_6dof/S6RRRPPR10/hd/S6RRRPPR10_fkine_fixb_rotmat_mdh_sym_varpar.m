% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2018-11-23 17:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRPPR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:40:10
% EndTime: 2018-11-23 17:40:10
% DurationCPUTime: 0.20s
% Computational Cost: add. (584->81), mult. (646->87), div. (0->0), fcn. (727->16), ass. (0->57)
t73 = pkin(4) + pkin(9);
t39 = sin(qJ(2));
t40 = sin(qJ(1));
t42 = cos(qJ(1));
t64 = pkin(6) - qJ(2);
t53 = cos(t64) / 0.2e1;
t63 = pkin(6) + qJ(2);
t57 = cos(t63);
t44 = t53 + t57 / 0.2e1;
t12 = t40 * t39 - t42 * t44;
t72 = t12 * pkin(9);
t14 = t42 * t39 + t40 * t44;
t71 = t14 * pkin(9);
t52 = sin(t63) / 0.2e1;
t56 = sin(t64);
t18 = t52 + t56 / 0.2e1;
t70 = t18 * pkin(9);
t36 = cos(pkin(11));
t69 = t36 * pkin(5) + t73;
t68 = cos(qJ(3));
t35 = sin(pkin(6));
t67 = t40 * t35;
t66 = t42 * t35;
t65 = cos(pkin(6));
t62 = pkin(7) + 0;
t61 = t35 * t68;
t60 = t65 * pkin(8) + t62;
t34 = sin(pkin(11));
t59 = pkin(5) * t34 + qJ(4);
t58 = t42 * pkin(1) + pkin(8) * t67 + 0;
t20 = t53 - t57 / 0.2e1;
t55 = t20 * pkin(2) + t60;
t19 = t52 - t56 / 0.2e1;
t41 = cos(qJ(2));
t15 = -t40 * t19 + t42 * t41;
t54 = t15 * pkin(2) + t58;
t38 = sin(qJ(3));
t11 = t20 * t68 + t65 * t38;
t51 = t11 * pkin(3) + t55;
t6 = t15 * t68 + t38 * t67;
t50 = t6 * pkin(3) + t54;
t49 = t40 * pkin(1) - pkin(8) * t66 + 0;
t13 = t42 * t19 + t40 * t41;
t48 = t13 * pkin(2) + t49;
t4 = t13 * t68 - t38 * t66;
t47 = t4 * pkin(3) + t48;
t5 = t15 * t38 - t40 * t61;
t46 = t5 * qJ(4) + t50;
t10 = t20 * t38 - t65 * t68;
t45 = t10 * qJ(4) + t51;
t3 = t13 * t38 + t42 * t61;
t43 = t3 * qJ(4) + t47;
t37 = -pkin(10) - qJ(5);
t33 = pkin(11) + qJ(6);
t29 = cos(t33);
t28 = sin(t33);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t42, -t40, 0, 0; t40, t42, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t15, -t14, t67, t58; t13, -t12, -t66, t49; t20, t18, t65, t60; 0, 0, 0, 1; t6, -t5, t14, t54 + t71; t4, -t3, t12, t48 + t72; t11, -t10, -t18, t55 - t70; 0, 0, 0, 1; t14, -t6, t5, t46 + t71; t12, -t4, t3, t43 + t72; -t18, -t11, t10, t45 - t70; 0, 0, 0, 1; t14 * t36 + t5 * t34, -t14 * t34 + t5 * t36, t6, t6 * qJ(5) + t73 * t14 + t46; t12 * t36 + t3 * t34, -t12 * t34 + t3 * t36, t4, t4 * qJ(5) + t73 * t12 + t43; t10 * t34 - t18 * t36, t10 * t36 + t18 * t34, t11, t11 * qJ(5) - t73 * t18 + t45; 0, 0, 0, 1; t14 * t29 + t5 * t28, -t14 * t28 + t5 * t29, t6, t69 * t14 - t6 * t37 + t59 * t5 + t50; t12 * t29 + t3 * t28, -t12 * t28 + t3 * t29, t4, t69 * t12 + t59 * t3 - t4 * t37 + t47; t10 * t28 - t18 * t29, t10 * t29 + t18 * t28, t11, t59 * t10 - t11 * t37 - t69 * t18 + t51; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
