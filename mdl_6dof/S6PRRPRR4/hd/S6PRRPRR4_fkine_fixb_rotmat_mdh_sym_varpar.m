% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:16:33
% EndTime: 2018-11-23 15:16:34
% DurationCPUTime: 0.20s
% Computational Cost: add. (671->76), mult. (789->82), div. (0->0), fcn. (900->16), ass. (0->59)
t79 = pkin(8) - pkin(9);
t41 = sin(pkin(11));
t43 = cos(pkin(11));
t47 = sin(qJ(2));
t71 = pkin(6) - qJ(2);
t61 = cos(t71) / 0.2e1;
t70 = pkin(6) + qJ(2);
t65 = cos(t70);
t55 = t61 + t65 / 0.2e1;
t22 = t41 * t47 - t43 * t55;
t78 = t22 * pkin(8);
t24 = t41 * t55 + t43 * t47;
t77 = t24 * pkin(8);
t60 = sin(t70) / 0.2e1;
t64 = sin(t71);
t30 = t60 + t64 / 0.2e1;
t76 = t30 * pkin(8);
t75 = cos(qJ(3));
t42 = sin(pkin(6));
t74 = t41 * t42;
t73 = t43 * t42;
t72 = cos(pkin(6));
t69 = qJ(1) + 0;
t68 = t42 * t75;
t67 = t43 * pkin(1) + pkin(7) * t74 + 0;
t66 = t72 * pkin(7) + t69;
t31 = t60 - t64 / 0.2e1;
t50 = cos(qJ(2));
t25 = -t41 * t31 + t43 * t50;
t63 = t25 * pkin(2) + t67;
t32 = t61 - t65 / 0.2e1;
t62 = t32 * pkin(2) + t66;
t59 = t41 * pkin(1) - pkin(7) * t73 + 0;
t23 = t43 * t31 + t41 * t50;
t58 = t23 * pkin(2) + t59;
t46 = sin(qJ(3));
t15 = t25 * t46 - t41 * t68;
t16 = t25 * t75 + t46 * t74;
t57 = t16 * pkin(3) + t15 * qJ(4) + t63;
t26 = t32 * t46 - t72 * t75;
t27 = t32 * t75 + t72 * t46;
t56 = t27 * pkin(3) + t26 * qJ(4) + t62;
t13 = t23 * t46 + t43 * t68;
t14 = t23 * t75 - t46 * t73;
t54 = t14 * pkin(3) + t13 * qJ(4) + t58;
t53 = t16 * pkin(4) + t79 * t24 + t57;
t52 = t27 * pkin(4) - t79 * t30 + t56;
t51 = t14 * pkin(4) + t79 * t22 + t54;
t49 = cos(qJ(5));
t48 = cos(qJ(6));
t45 = sin(qJ(5));
t44 = sin(qJ(6));
t6 = t26 * t45 + t27 * t49;
t5 = -t26 * t49 + t27 * t45;
t4 = t15 * t45 + t16 * t49;
t3 = -t15 * t49 + t16 * t45;
t2 = t13 * t45 + t14 * t49;
t1 = -t13 * t49 + t14 * t45;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t43, -t41, 0, 0; t41, t43, 0, 0; 0, 0, 1, t69; 0, 0, 0, 1; t25, -t24, t74, t67; t23, -t22, -t73, t59; t32, t30, t72, t66; 0, 0, 0, 1; t16, -t15, t24, t63 + t77; t14, -t13, t22, t58 + t78; t27, -t26, -t30, t62 - t76; 0, 0, 0, 1; t16, t24, t15, t57 + t77; t14, t22, t13, t54 + t78; t27, -t30, t26, t56 - t76; 0, 0, 0, 1; t4, -t3, -t24, t53; t2, -t1, -t22, t51; t6, -t5, t30, t52; 0, 0, 0, 1; -t24 * t44 + t4 * t48, -t24 * t48 - t4 * t44, t3, t4 * pkin(5) + t3 * pkin(10) + t53; t2 * t48 - t22 * t44, -t2 * t44 - t22 * t48, t1, t2 * pkin(5) + t1 * pkin(10) + t51; t30 * t44 + t6 * t48, t30 * t48 - t6 * t44, t5, t6 * pkin(5) + t5 * pkin(10) + t52; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
