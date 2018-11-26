% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRPRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:44:20
% EndTime: 2018-11-23 17:44:20
% DurationCPUTime: 0.18s
% Computational Cost: add. (572->79), mult. (561->90), div. (0->0), fcn. (633->16), ass. (0->59)
t48 = sin(qJ(2));
t49 = sin(qJ(1));
t53 = cos(qJ(1));
t69 = pkin(6) - qJ(2);
t60 = cos(t69) / 0.2e1;
t68 = pkin(6) + qJ(2);
t62 = cos(t68);
t55 = t60 + t62 / 0.2e1;
t17 = t49 * t48 - t53 * t55;
t46 = sin(qJ(5));
t75 = t17 * t46;
t19 = t53 * t48 + t49 * t55;
t74 = t19 * t46;
t59 = sin(t68) / 0.2e1;
t61 = sin(t69);
t24 = t59 + t61 / 0.2e1;
t73 = t24 * t46;
t43 = cos(pkin(6));
t47 = sin(qJ(3));
t72 = t43 * t47;
t42 = sin(pkin(6));
t71 = t49 * t42;
t70 = t53 * t42;
t67 = pkin(7) + 0;
t66 = t49 * pkin(1) + 0;
t65 = t47 * t71;
t64 = t43 * pkin(8) + t67;
t63 = t53 * pkin(1) + pkin(8) * t71 + 0;
t58 = -pkin(8) * t70 + t66;
t26 = t60 - t62 / 0.2e1;
t51 = cos(qJ(3));
t35 = t51 * pkin(3) + pkin(2);
t45 = -qJ(4) - pkin(9);
t57 = pkin(3) * t72 + t24 * t45 + t26 * t35 + t64;
t25 = t59 - t61 / 0.2e1;
t52 = cos(qJ(2));
t20 = -t49 * t25 + t53 * t52;
t56 = pkin(3) * t65 - t19 * t45 + t20 * t35 + t63;
t18 = t53 * t25 + t49 * t52;
t54 = t18 * t35 - t17 * t45 + (-pkin(3) * t47 - pkin(8)) * t70 + t66;
t50 = cos(qJ(5));
t44 = -qJ(6) - pkin(10);
t41 = qJ(3) + pkin(11);
t37 = cos(t41);
t36 = sin(t41);
t34 = t50 * pkin(5) + pkin(4);
t14 = t26 * t37 + t43 * t36;
t13 = t26 * t36 - t43 * t37;
t10 = t20 * t37 + t36 * t71;
t9 = t20 * t36 - t37 * t71;
t8 = t18 * t37 - t36 * t70;
t7 = t18 * t36 + t37 * t70;
t6 = t14 * t50 - t73;
t5 = -t14 * t46 - t24 * t50;
t4 = t10 * t50 + t74;
t3 = -t10 * t46 + t19 * t50;
t2 = t8 * t50 + t75;
t1 = t17 * t50 - t8 * t46;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t53, -t49, 0, 0; t49, t53, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t20, -t19, t71, t63; t18, -t17, -t70, t58; t26, t24, t43, t64; 0, 0, 0, 1; t20 * t51 + t65, -t20 * t47 + t51 * t71, t19, t20 * pkin(2) + t19 * pkin(9) + t63; t18 * t51 - t47 * t70, -t18 * t47 - t51 * t70, t17, t18 * pkin(2) + t17 * pkin(9) + t58; t26 * t51 + t72, -t26 * t47 + t43 * t51, -t24, t26 * pkin(2) - t24 * pkin(9) + t64; 0, 0, 0, 1; t10, -t9, t19, t56; t8, -t7, t17, t54; t14, -t13, -t24, t57; 0, 0, 0, 1; t4, t3, t9, t10 * pkin(4) + t9 * pkin(10) + t56; t2, t1, t7, t8 * pkin(4) + t7 * pkin(10) + t54; t6, t5, t13, t14 * pkin(4) + t13 * pkin(10) + t57; 0, 0, 0, 1; t4, t3, t9, pkin(5) * t74 + t10 * t34 - t9 * t44 + t56; t2, t1, t7, pkin(5) * t75 + t8 * t34 - t7 * t44 + t54; t6, t5, t13, -pkin(5) * t73 - t13 * t44 + t14 * t34 + t57; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
