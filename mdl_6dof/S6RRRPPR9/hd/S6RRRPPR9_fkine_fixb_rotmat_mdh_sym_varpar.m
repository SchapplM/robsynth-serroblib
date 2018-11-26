% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2018-11-23 17:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRPPR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:39:12
% EndTime: 2018-11-23 17:39:13
% DurationCPUTime: 0.17s
% Computational Cost: add. (752->76), mult. (877->82), div. (0->0), fcn. (997->16), ass. (0->59)
t68 = pkin(6) + qJ(2);
t79 = sin(t68) / 0.2e1;
t78 = cos(qJ(3));
t42 = sin(pkin(6));
t46 = sin(qJ(1));
t77 = t46 * t42;
t49 = cos(qJ(1));
t76 = t49 * t42;
t75 = -pkin(10) + qJ(4);
t69 = pkin(6) - qJ(2);
t62 = sin(t69);
t31 = t79 - t62 / 0.2e1;
t48 = cos(qJ(2));
t24 = t49 * t31 + t46 * t48;
t44 = sin(qJ(3));
t66 = t42 * t78;
t12 = t24 * t44 + t49 * t66;
t74 = t12 * qJ(4);
t26 = -t46 * t31 + t49 * t48;
t14 = t26 * t44 - t46 * t66;
t73 = t14 * qJ(4);
t61 = cos(t69) / 0.2e1;
t63 = cos(t68);
t32 = t61 - t63 / 0.2e1;
t71 = cos(pkin(6));
t21 = t32 * t44 - t71 * t78;
t72 = t21 * qJ(4);
t70 = cos(pkin(11));
t67 = pkin(7) + 0;
t65 = t71 * pkin(8) + t67;
t64 = t49 * pkin(1) + pkin(8) * t77 + 0;
t60 = t46 * pkin(1) - pkin(8) * t76 + 0;
t30 = t79 + t62 / 0.2e1;
t59 = t32 * pkin(2) - t30 * pkin(9) + t65;
t45 = sin(qJ(2));
t54 = t61 + t63 / 0.2e1;
t25 = t49 * t45 + t46 * t54;
t58 = t26 * pkin(2) + t25 * pkin(9) + t64;
t22 = t32 * t78 + t71 * t44;
t57 = t22 * pkin(3) + t59;
t15 = t26 * t78 + t44 * t77;
t56 = t15 * pkin(3) + t58;
t23 = t46 * t45 - t49 * t54;
t55 = t24 * pkin(2) + t23 * pkin(9) + t60;
t13 = t24 * t78 - t44 * t76;
t53 = t13 * pkin(3) + t55;
t41 = sin(pkin(11));
t8 = t22 * t41 + t30 * t70;
t9 = t22 * t70 - t30 * t41;
t52 = t9 * pkin(4) + t8 * qJ(5) + t57;
t5 = t15 * t41 - t25 * t70;
t6 = t15 * t70 + t25 * t41;
t51 = t6 * pkin(4) + t5 * qJ(5) + t56;
t3 = t13 * t41 - t23 * t70;
t4 = t13 * t70 + t23 * t41;
t50 = t4 * pkin(4) + t3 * qJ(5) + t53;
t47 = cos(qJ(6));
t43 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t49, -t46, 0, 0; t46, t49, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t26, -t25, t77, t64; t24, -t23, -t76, t60; t32, t30, t71, t65; 0, 0, 0, 1; t15, -t14, t25, t58; t13, -t12, t23, t55; t22, -t21, -t30, t59; 0, 0, 0, 1; t6, -t5, t14, t56 + t73; t4, -t3, t12, t53 + t74; t9, -t8, t21, t57 + t72; 0, 0, 0, 1; t6, t14, t5, t51 + t73; t4, t12, t3, t50 + t74; t9, t21, t8, t52 + t72; 0, 0, 0, 1; t5 * t43 + t6 * t47, -t6 * t43 + t5 * t47, -t14, t6 * pkin(5) + t75 * t14 + t51; t3 * t43 + t4 * t47, t3 * t47 - t4 * t43, -t12, t4 * pkin(5) + t75 * t12 + t50; t8 * t43 + t9 * t47, -t9 * t43 + t8 * t47, -t21, t9 * pkin(5) + t75 * t21 + t52; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
