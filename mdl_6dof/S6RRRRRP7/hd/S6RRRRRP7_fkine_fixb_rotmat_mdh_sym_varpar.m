% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP7
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
% Datum: 2018-11-23 18:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRRRP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:30:44
% EndTime: 2018-11-23 18:30:44
% DurationCPUTime: 0.19s
% Computational Cost: add. (572->79), mult. (561->90), div. (0->0), fcn. (633->16), ass. (0->59)
t47 = sin(qJ(2));
t48 = sin(qJ(1));
t52 = cos(qJ(1));
t69 = pkin(6) - qJ(2);
t60 = cos(t69) / 0.2e1;
t68 = pkin(6) + qJ(2);
t62 = cos(t68);
t55 = t60 + t62 / 0.2e1;
t17 = t48 * t47 - t52 * t55;
t45 = sin(qJ(5));
t75 = t17 * t45;
t19 = t52 * t47 + t48 * t55;
t74 = t19 * t45;
t59 = sin(t68) / 0.2e1;
t61 = sin(t69);
t24 = t59 + t61 / 0.2e1;
t73 = t24 * t45;
t43 = cos(pkin(6));
t46 = sin(qJ(3));
t72 = t43 * t46;
t42 = sin(pkin(6));
t71 = t48 * t42;
t70 = t52 * t42;
t67 = pkin(7) + 0;
t66 = t48 * pkin(1) + 0;
t65 = t46 * t71;
t64 = t43 * pkin(8) + t67;
t63 = t52 * pkin(1) + pkin(8) * t71 + 0;
t58 = -pkin(8) * t70 + t66;
t26 = t60 - t62 / 0.2e1;
t50 = cos(qJ(3));
t35 = t50 * pkin(3) + pkin(2);
t53 = -pkin(10) - pkin(9);
t57 = pkin(3) * t72 + t24 * t53 + t26 * t35 + t64;
t25 = t59 - t61 / 0.2e1;
t51 = cos(qJ(2));
t20 = -t48 * t25 + t52 * t51;
t56 = pkin(3) * t65 - t19 * t53 + t20 * t35 + t63;
t18 = t52 * t25 + t48 * t51;
t54 = t18 * t35 - t17 * t53 + (-pkin(3) * t46 - pkin(8)) * t70 + t66;
t49 = cos(qJ(5));
t44 = -qJ(6) - pkin(11);
t41 = qJ(3) + qJ(4);
t37 = cos(t41);
t36 = sin(t41);
t34 = t49 * pkin(5) + pkin(4);
t14 = t26 * t37 + t43 * t36;
t13 = t26 * t36 - t43 * t37;
t10 = t20 * t37 + t36 * t71;
t9 = t20 * t36 - t37 * t71;
t8 = t18 * t37 - t36 * t70;
t7 = t18 * t36 + t37 * t70;
t6 = t14 * t49 - t73;
t5 = -t14 * t45 - t24 * t49;
t4 = t10 * t49 + t74;
t3 = -t10 * t45 + t19 * t49;
t2 = t8 * t49 + t75;
t1 = t17 * t49 - t8 * t45;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t52, -t48, 0, 0; t48, t52, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t20, -t19, t71, t63; t18, -t17, -t70, t58; t26, t24, t43, t64; 0, 0, 0, 1; t20 * t50 + t65, -t20 * t46 + t50 * t71, t19, t20 * pkin(2) + t19 * pkin(9) + t63; t18 * t50 - t46 * t70, -t18 * t46 - t50 * t70, t17, t18 * pkin(2) + t17 * pkin(9) + t58; t26 * t50 + t72, -t26 * t46 + t43 * t50, -t24, t26 * pkin(2) - t24 * pkin(9) + t64; 0, 0, 0, 1; t10, -t9, t19, t56; t8, -t7, t17, t54; t14, -t13, -t24, t57; 0, 0, 0, 1; t4, t3, t9, t10 * pkin(4) + t9 * pkin(11) + t56; t2, t1, t7, t8 * pkin(4) + t7 * pkin(11) + t54; t6, t5, t13, t14 * pkin(4) + t13 * pkin(11) + t57; 0, 0, 0, 1; t4, t3, t9, pkin(5) * t74 + t10 * t34 - t9 * t44 + t56; t2, t1, t7, pkin(5) * t75 + t8 * t34 - t7 * t44 + t54; t6, t5, t13, -pkin(5) * t73 - t13 * t44 + t14 * t34 + t57; 0, 0, 0, 1;];
T_ges = t11;
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
