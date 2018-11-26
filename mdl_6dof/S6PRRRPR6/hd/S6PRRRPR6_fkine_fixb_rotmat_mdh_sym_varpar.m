% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2018-11-23 15:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRRPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:25:17
% EndTime: 2018-11-23 15:25:17
% DurationCPUTime: 0.20s
% Computational Cost: add. (752->76), mult. (877->82), div. (0->0), fcn. (997->16), ass. (0->59)
t68 = pkin(6) + qJ(2);
t79 = sin(t68) / 0.2e1;
t78 = pkin(9) - pkin(10);
t69 = pkin(6) - qJ(2);
t62 = sin(t69);
t31 = t79 - t62 / 0.2e1;
t41 = sin(pkin(11));
t43 = cos(pkin(11));
t49 = cos(qJ(2));
t22 = t43 * t31 + t41 * t49;
t46 = sin(qJ(3));
t42 = sin(pkin(6));
t74 = cos(qJ(3));
t66 = t42 * t74;
t12 = t22 * t46 + t43 * t66;
t77 = t12 * pkin(9);
t24 = -t41 * t31 + t43 * t49;
t14 = t24 * t46 - t41 * t66;
t76 = t14 * pkin(9);
t61 = cos(t69) / 0.2e1;
t63 = cos(t68);
t32 = t61 - t63 / 0.2e1;
t70 = cos(pkin(6));
t25 = t32 * t46 - t70 * t74;
t75 = t25 * pkin(9);
t73 = cos(qJ(4));
t72 = t41 * t42;
t71 = t43 * t42;
t67 = qJ(1) + 0;
t65 = t43 * pkin(1) + pkin(7) * t72 + 0;
t64 = t70 * pkin(7) + t67;
t60 = t41 * pkin(1) - pkin(7) * t71 + 0;
t47 = sin(qJ(2));
t54 = t61 + t63 / 0.2e1;
t23 = t41 * t54 + t43 * t47;
t59 = t24 * pkin(2) + t23 * pkin(8) + t65;
t30 = t79 + t62 / 0.2e1;
t58 = t32 * pkin(2) - t30 * pkin(8) + t64;
t15 = t24 * t74 + t46 * t72;
t57 = t15 * pkin(3) + t59;
t26 = t32 * t74 + t70 * t46;
t56 = t26 * pkin(3) + t58;
t21 = t41 * t47 - t43 * t54;
t55 = t22 * pkin(2) + t21 * pkin(8) + t60;
t13 = t22 * t74 - t46 * t71;
t53 = t13 * pkin(3) + t55;
t45 = sin(qJ(4));
t5 = t15 * t45 - t23 * t73;
t6 = t15 * t73 + t23 * t45;
t52 = t6 * pkin(4) + t5 * qJ(5) + t57;
t8 = t26 * t45 + t30 * t73;
t9 = t26 * t73 - t30 * t45;
t51 = t9 * pkin(4) + t8 * qJ(5) + t56;
t3 = t13 * t45 - t21 * t73;
t4 = t13 * t73 + t21 * t45;
t50 = t4 * pkin(4) + t3 * qJ(5) + t53;
t48 = cos(qJ(6));
t44 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t43, -t41, 0, 0; t41, t43, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t24, -t23, t72, t65; t22, -t21, -t71, t60; t32, t30, t70, t64; 0, 0, 0, 1; t15, -t14, t23, t59; t13, -t12, t21, t55; t26, -t25, -t30, t58; 0, 0, 0, 1; t6, -t5, t14, t57 + t76; t4, -t3, t12, t53 + t77; t9, -t8, t25, t56 + t75; 0, 0, 0, 1; t6, t14, t5, t52 + t76; t4, t12, t3, t50 + t77; t9, t25, t8, t51 + t75; 0, 0, 0, 1; t5 * t44 + t6 * t48, -t6 * t44 + t5 * t48, -t14, t6 * pkin(5) + t78 * t14 + t52; t3 * t44 + t4 * t48, t3 * t48 - t4 * t44, -t12, t4 * pkin(5) + t78 * t12 + t50; t8 * t44 + t9 * t48, -t9 * t44 + t8 * t48, -t25, t9 * pkin(5) + t78 * t25 + t51; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
