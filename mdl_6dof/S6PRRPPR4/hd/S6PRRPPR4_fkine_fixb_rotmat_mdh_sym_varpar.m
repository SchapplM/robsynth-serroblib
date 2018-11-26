% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRPPR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:09:58
% EndTime: 2018-11-23 15:09:58
% DurationCPUTime: 0.20s
% Computational Cost: add. (752->76), mult. (877->82), div. (0->0), fcn. (997->16), ass. (0->59)
t68 = pkin(6) + qJ(2);
t79 = sin(t68) / 0.2e1;
t78 = cos(qJ(3));
t42 = sin(pkin(10));
t43 = sin(pkin(6));
t77 = t42 * t43;
t44 = cos(pkin(10));
t76 = t44 * t43;
t75 = -pkin(9) + qJ(4);
t69 = pkin(6) - qJ(2);
t62 = sin(t69);
t31 = t79 - t62 / 0.2e1;
t49 = cos(qJ(2));
t22 = t44 * t31 + t42 * t49;
t46 = sin(qJ(3));
t66 = t43 * t78;
t12 = t22 * t46 + t44 * t66;
t74 = t12 * qJ(4);
t24 = -t42 * t31 + t44 * t49;
t14 = t24 * t46 - t42 * t66;
t73 = t14 * qJ(4);
t61 = cos(t69) / 0.2e1;
t63 = cos(t68);
t32 = t61 - t63 / 0.2e1;
t71 = cos(pkin(6));
t25 = t32 * t46 - t71 * t78;
t72 = t25 * qJ(4);
t70 = cos(pkin(11));
t67 = qJ(1) + 0;
t65 = t44 * pkin(1) + pkin(7) * t77 + 0;
t64 = t71 * pkin(7) + t67;
t60 = t42 * pkin(1) - pkin(7) * t76 + 0;
t47 = sin(qJ(2));
t54 = t61 + t63 / 0.2e1;
t23 = t42 * t54 + t44 * t47;
t59 = t24 * pkin(2) + t23 * pkin(8) + t65;
t30 = t79 + t62 / 0.2e1;
t58 = t32 * pkin(2) - t30 * pkin(8) + t64;
t15 = t24 * t78 + t46 * t77;
t57 = t15 * pkin(3) + t59;
t26 = t32 * t78 + t71 * t46;
t56 = t26 * pkin(3) + t58;
t21 = t42 * t47 - t44 * t54;
t55 = t22 * pkin(2) + t21 * pkin(8) + t60;
t13 = t22 * t78 - t46 * t76;
t53 = t13 * pkin(3) + t55;
t41 = sin(pkin(11));
t5 = t15 * t41 - t23 * t70;
t6 = t15 * t70 + t23 * t41;
t52 = t6 * pkin(4) + t5 * qJ(5) + t57;
t8 = t26 * t41 + t30 * t70;
t9 = t26 * t70 - t30 * t41;
t51 = t9 * pkin(4) + t8 * qJ(5) + t56;
t3 = t13 * t41 - t21 * t70;
t4 = t13 * t70 + t21 * t41;
t50 = t4 * pkin(4) + t3 * qJ(5) + t53;
t48 = cos(qJ(6));
t45 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t44, -t42, 0, 0; t42, t44, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t24, -t23, t77, t65; t22, -t21, -t76, t60; t32, t30, t71, t64; 0, 0, 0, 1; t15, -t14, t23, t59; t13, -t12, t21, t55; t26, -t25, -t30, t58; 0, 0, 0, 1; t6, -t5, t14, t57 + t73; t4, -t3, t12, t53 + t74; t9, -t8, t25, t56 + t72; 0, 0, 0, 1; t6, t14, t5, t52 + t73; t4, t12, t3, t50 + t74; t9, t25, t8, t51 + t72; 0, 0, 0, 1; t5 * t45 + t6 * t48, -t6 * t45 + t5 * t48, -t14, t6 * pkin(5) + t75 * t14 + t52; t3 * t45 + t4 * t48, t3 * t48 - t4 * t45, -t12, t4 * pkin(5) + t75 * t12 + t50; t8 * t45 + t9 * t48, -t9 * t45 + t8 * t48, -t25, t9 * pkin(5) + t75 * t25 + t51; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
