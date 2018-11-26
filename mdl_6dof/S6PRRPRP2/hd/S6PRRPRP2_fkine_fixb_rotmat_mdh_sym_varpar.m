% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2018-11-23 15:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRPRP2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:11:32
% EndTime: 2018-11-23 15:11:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (626->77), mult. (616->86), div. (0->0), fcn. (697->16), ass. (0->57)
t47 = sin(pkin(10));
t48 = sin(pkin(6));
t79 = t47 * t48;
t49 = cos(pkin(10));
t78 = t49 * t48;
t50 = cos(pkin(6));
t53 = sin(qJ(3));
t77 = t50 * t53;
t76 = pkin(6) - qJ(2);
t75 = pkin(6) + qJ(2);
t74 = t47 * pkin(1) + 0;
t73 = t53 * t79;
t72 = qJ(1) + 0;
t71 = t49 * pkin(1) + pkin(7) * t79 + 0;
t70 = t50 * pkin(7) + t72;
t69 = cos(t75);
t68 = sin(t76);
t67 = cos(t76) / 0.2e1;
t66 = sin(t75) / 0.2e1;
t65 = -pkin(7) * t78 + t74;
t54 = sin(qJ(2));
t60 = t67 + t69 / 0.2e1;
t24 = t47 * t60 + t49 * t54;
t31 = t66 - t68 / 0.2e1;
t57 = cos(qJ(2));
t25 = -t47 * t31 + t49 * t57;
t56 = cos(qJ(3));
t40 = t56 * pkin(3) + pkin(2);
t51 = -qJ(4) - pkin(8);
t64 = pkin(3) * t73 - t24 * t51 + t25 * t40 + t71;
t30 = t66 + t68 / 0.2e1;
t32 = t67 - t69 / 0.2e1;
t63 = pkin(3) * t77 + t30 * t51 + t32 * t40 + t70;
t46 = qJ(3) + pkin(11);
t41 = sin(t46);
t42 = cos(t46);
t11 = t25 * t41 - t42 * t79;
t12 = t25 * t42 + t41 * t79;
t62 = t12 * pkin(4) + t11 * pkin(9) + t64;
t16 = t32 * t41 - t50 * t42;
t17 = t32 * t42 + t50 * t41;
t61 = t17 * pkin(4) + t16 * pkin(9) + t63;
t22 = t47 * t54 - t49 * t60;
t23 = t49 * t31 + t47 * t57;
t59 = t23 * t40 - t22 * t51 + (-pkin(3) * t53 - pkin(7)) * t78 + t74;
t10 = t23 * t42 - t41 * t78;
t9 = t23 * t41 + t42 * t78;
t58 = t10 * pkin(4) + t9 * pkin(9) + t59;
t55 = cos(qJ(5));
t52 = sin(qJ(5));
t6 = t17 * t55 - t30 * t52;
t5 = t17 * t52 + t30 * t55;
t4 = t12 * t55 + t24 * t52;
t3 = t12 * t52 - t24 * t55;
t2 = t10 * t55 + t22 * t52;
t1 = t10 * t52 - t22 * t55;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t49, -t47, 0, 0; t47, t49, 0, 0; 0, 0, 1, t72; 0, 0, 0, 1; t25, -t24, t79, t71; t23, -t22, -t78, t65; t32, t30, t50, t70; 0, 0, 0, 1; t25 * t56 + t73, -t25 * t53 + t56 * t79, t24, t25 * pkin(2) + t24 * pkin(8) + t71; t23 * t56 - t53 * t78, -t23 * t53 - t56 * t78, t22, t23 * pkin(2) + t22 * pkin(8) + t65; t32 * t56 + t77, -t32 * t53 + t50 * t56, -t30, t32 * pkin(2) - t30 * pkin(8) + t70; 0, 0, 0, 1; t12, -t11, t24, t64; t10, -t9, t22, t59; t17, -t16, -t30, t63; 0, 0, 0, 1; t4, -t3, t11, t62; t2, -t1, t9, t58; t6, -t5, t16, t61; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(5) + t3 * qJ(6) + t62; t2, t9, t1, t2 * pkin(5) + t1 * qJ(6) + t58; t6, t16, t5, t6 * pkin(5) + t5 * qJ(6) + t61; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
