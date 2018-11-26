% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:42
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:41:23
% EndTime: 2018-11-23 18:41:23
% DurationCPUTime: 0.19s
% Computational Cost: add. (584->86), mult. (561->102), div. (0->0), fcn. (633->18), ass. (0->56)
t43 = sin(qJ(2));
t44 = sin(qJ(1));
t48 = cos(qJ(1));
t66 = pkin(6) - qJ(2);
t57 = cos(t66) / 0.2e1;
t65 = pkin(6) + qJ(2);
t59 = cos(t65);
t52 = t57 + t59 / 0.2e1;
t11 = t44 * t43 - t48 * t52;
t41 = sin(qJ(5));
t72 = t11 * t41;
t13 = t48 * t43 + t44 * t52;
t71 = t13 * t41;
t56 = sin(t65) / 0.2e1;
t58 = sin(t66);
t18 = t56 + t58 / 0.2e1;
t70 = t18 * t41;
t40 = cos(pkin(6));
t42 = sin(qJ(3));
t69 = t40 * t42;
t39 = sin(pkin(6));
t68 = t44 * t39;
t67 = t48 * t39;
t64 = pkin(7) + 0;
t63 = t44 * pkin(1) + 0;
t62 = t42 * t68;
t61 = t40 * pkin(8) + t64;
t60 = t48 * pkin(1) + pkin(8) * t68 + 0;
t55 = -pkin(8) * t67 + t63;
t20 = t57 - t59 / 0.2e1;
t46 = cos(qJ(3));
t29 = t46 * pkin(3) + pkin(2);
t50 = -pkin(10) - pkin(9);
t54 = pkin(3) * t69 + t18 * t50 + t20 * t29 + t61;
t19 = t56 - t58 / 0.2e1;
t47 = cos(qJ(2));
t14 = -t44 * t19 + t48 * t47;
t53 = pkin(3) * t62 - t13 * t50 + t14 * t29 + t60;
t12 = t48 * t19 + t44 * t47;
t51 = t12 * t29 - t11 * t50 + (-pkin(3) * t42 - pkin(8)) * t67 + t63;
t49 = -pkin(12) - pkin(11);
t45 = cos(qJ(5));
t38 = qJ(3) + qJ(4);
t37 = qJ(5) + qJ(6);
t33 = cos(t38);
t32 = cos(t37);
t31 = sin(t38);
t30 = sin(t37);
t28 = t45 * pkin(5) + pkin(4);
t8 = t20 * t33 + t40 * t31;
t7 = t20 * t31 - t40 * t33;
t4 = t14 * t33 + t31 * t68;
t3 = t14 * t31 - t33 * t68;
t2 = t12 * t33 - t31 * t67;
t1 = t12 * t31 + t33 * t67;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t48, -t44, 0, 0; t44, t48, 0, 0; 0, 0, 1, t64; 0, 0, 0, 1; t14, -t13, t68, t60; t12, -t11, -t67, t55; t20, t18, t40, t61; 0, 0, 0, 1; t14 * t46 + t62, -t14 * t42 + t46 * t68, t13, t14 * pkin(2) + t13 * pkin(9) + t60; t12 * t46 - t42 * t67, -t12 * t42 - t46 * t67, t11, t12 * pkin(2) + t11 * pkin(9) + t55; t20 * t46 + t69, -t20 * t42 + t40 * t46, -t18, t20 * pkin(2) - t18 * pkin(9) + t61; 0, 0, 0, 1; t4, -t3, t13, t53; t2, -t1, t11, t51; t8, -t7, -t18, t54; 0, 0, 0, 1; t4 * t45 + t71, t13 * t45 - t4 * t41, t3, t4 * pkin(4) + t3 * pkin(11) + t53; t2 * t45 + t72, t11 * t45 - t2 * t41, t1, t2 * pkin(4) + t1 * pkin(11) + t51; t8 * t45 - t70, -t18 * t45 - t8 * t41, t7, t8 * pkin(4) + t7 * pkin(11) + t54; 0, 0, 0, 1; t13 * t30 + t4 * t32, t13 * t32 - t4 * t30, t3, pkin(5) * t71 + t4 * t28 - t3 * t49 + t53; t11 * t30 + t2 * t32, t11 * t32 - t2 * t30, t1, pkin(5) * t72 - t1 * t49 + t2 * t28 + t51; -t18 * t30 + t8 * t32, -t18 * t32 - t8 * t30, t7, -pkin(5) * t70 + t8 * t28 - t7 * t49 + t54; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
