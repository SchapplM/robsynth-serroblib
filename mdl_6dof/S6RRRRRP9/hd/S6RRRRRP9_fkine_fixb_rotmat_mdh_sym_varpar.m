% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRP9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:32:44
% EndTime: 2018-11-23 18:32:44
% DurationCPUTime: 0.21s
% Computational Cost: add. (597->79), mult. (635->92), div. (0->0), fcn. (717->16), ass. (0->59)
t41 = sin(qJ(4));
t64 = t41 * pkin(4) + pkin(9);
t48 = -pkin(11) - pkin(10);
t39 = qJ(4) + qJ(5);
t32 = sin(t39);
t72 = pkin(5) * t32 + t64;
t45 = cos(qJ(4));
t31 = t45 * pkin(4) + pkin(3);
t71 = cos(qJ(3));
t40 = sin(pkin(6));
t44 = sin(qJ(1));
t70 = t44 * t40;
t47 = cos(qJ(1));
t69 = t47 * t40;
t68 = cos(pkin(6));
t67 = pkin(6) - qJ(2);
t66 = pkin(6) + qJ(2);
t65 = pkin(7) + 0;
t63 = t40 * t71;
t62 = t68 * pkin(8) + t65;
t61 = t47 * pkin(1) + pkin(8) * t70 + 0;
t60 = cos(t66);
t59 = sin(t67);
t56 = cos(t67) / 0.2e1;
t24 = t56 - t60 / 0.2e1;
t58 = t24 * pkin(2) + t62;
t55 = sin(t66) / 0.2e1;
t23 = t55 - t59 / 0.2e1;
t46 = cos(qJ(2));
t18 = -t44 * t23 + t47 * t46;
t57 = t18 * pkin(2) + t61;
t54 = t44 * pkin(1) - pkin(8) * t69 + 0;
t22 = t55 + t59 / 0.2e1;
t53 = -t22 * pkin(9) + t58;
t16 = t47 * t23 + t44 * t46;
t52 = t16 * pkin(2) + t54;
t43 = sin(qJ(2));
t49 = t56 + t60 / 0.2e1;
t17 = t47 * t43 + t44 * t49;
t51 = t17 * pkin(9) + t57;
t15 = t44 * t43 - t47 * t49;
t50 = t15 * pkin(9) + t52;
t42 = sin(qJ(3));
t38 = -qJ(6) + t48;
t33 = cos(t39);
t21 = pkin(5) * t33 + t31;
t14 = t24 * t71 + t68 * t42;
t13 = t24 * t42 - t68 * t71;
t10 = t18 * t71 + t42 * t70;
t9 = t18 * t42 - t44 * t63;
t8 = t16 * t71 - t42 * t69;
t7 = t16 * t42 + t47 * t63;
t6 = t14 * t33 - t22 * t32;
t5 = -t14 * t32 - t22 * t33;
t4 = t10 * t33 + t17 * t32;
t3 = -t10 * t32 + t17 * t33;
t2 = t15 * t32 + t8 * t33;
t1 = t15 * t33 - t8 * t32;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t47, -t44, 0, 0; t44, t47, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t18, -t17, t70, t61; t16, -t15, -t69, t54; t24, t22, t68, t62; 0, 0, 0, 1; t10, -t9, t17, t51; t8, -t7, t15, t50; t14, -t13, -t22, t53; 0, 0, 0, 1; t10 * t45 + t17 * t41, -t10 * t41 + t17 * t45, t9, t10 * pkin(3) + t9 * pkin(10) + t51; t15 * t41 + t8 * t45, t15 * t45 - t8 * t41, t7, t8 * pkin(3) + t7 * pkin(10) + t50; t14 * t45 - t22 * t41, -t14 * t41 - t22 * t45, t13, t14 * pkin(3) + t13 * pkin(10) + t53; 0, 0, 0, 1; t4, t3, t9, t10 * t31 + t64 * t17 - t9 * t48 + t57; t2, t1, t7, t64 * t15 + t8 * t31 - t7 * t48 + t52; t6, t5, t13, -t13 * t48 + t14 * t31 - t64 * t22 + t58; 0, 0, 0, 1; t4, t3, t9, t10 * t21 + t72 * t17 - t9 * t38 + t57; t2, t1, t7, t72 * t15 + t8 * t21 - t7 * t38 + t52; t6, t5, t13, -t13 * t38 + t14 * t21 - t72 * t22 + t58; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
