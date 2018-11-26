% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRP10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:17:31
% EndTime: 2018-11-23 17:17:31
% DurationCPUTime: 0.21s
% Computational Cost: add. (626->77), mult. (616->86), div. (0->0), fcn. (697->16), ass. (0->57)
t47 = sin(pkin(11));
t50 = cos(pkin(6));
t79 = t50 * t47;
t48 = sin(pkin(6));
t54 = sin(qJ(1));
t78 = t54 * t48;
t57 = cos(qJ(1));
t77 = t57 * t48;
t76 = pkin(6) - qJ(2);
t75 = pkin(6) + qJ(2);
t74 = pkin(7) + 0;
t73 = t54 * pkin(1) + 0;
t72 = t47 * t78;
t71 = t50 * pkin(8) + t74;
t70 = t57 * pkin(1) + pkin(8) * t78 + 0;
t69 = cos(t75);
t68 = sin(t76);
t67 = cos(t76) / 0.2e1;
t66 = sin(t75) / 0.2e1;
t65 = -pkin(8) * t77 + t73;
t30 = t66 + t68 / 0.2e1;
t32 = t67 - t69 / 0.2e1;
t49 = cos(pkin(11));
t40 = t49 * pkin(3) + pkin(2);
t51 = -pkin(9) - qJ(3);
t64 = pkin(3) * t79 + t30 * t51 + t32 * t40 + t71;
t53 = sin(qJ(2));
t60 = t67 + t69 / 0.2e1;
t24 = t57 * t53 + t54 * t60;
t31 = t66 - t68 / 0.2e1;
t56 = cos(qJ(2));
t25 = -t54 * t31 + t57 * t56;
t63 = pkin(3) * t72 - t24 * t51 + t25 * t40 + t70;
t46 = pkin(11) + qJ(4);
t41 = sin(t46);
t42 = cos(t46);
t16 = t32 * t41 - t50 * t42;
t17 = t32 * t42 + t50 * t41;
t62 = t17 * pkin(4) + t16 * pkin(10) + t64;
t11 = t25 * t41 - t42 * t78;
t12 = t25 * t42 + t41 * t78;
t61 = t12 * pkin(4) + t11 * pkin(10) + t63;
t22 = t54 * t53 - t57 * t60;
t23 = t57 * t31 + t54 * t56;
t59 = t23 * t40 - t22 * t51 + (-pkin(3) * t47 - pkin(8)) * t77 + t73;
t10 = t23 * t42 - t41 * t77;
t9 = t23 * t41 + t42 * t77;
t58 = t10 * pkin(4) + t9 * pkin(10) + t59;
t55 = cos(qJ(5));
t52 = sin(qJ(5));
t6 = t17 * t55 - t30 * t52;
t5 = t17 * t52 + t30 * t55;
t4 = t12 * t55 + t24 * t52;
t3 = t12 * t52 - t24 * t55;
t2 = t10 * t55 + t22 * t52;
t1 = t10 * t52 - t22 * t55;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t57, -t54, 0, 0; t54, t57, 0, 0; 0, 0, 1, t74; 0, 0, 0, 1; t25, -t24, t78, t70; t23, -t22, -t77, t65; t32, t30, t50, t71; 0, 0, 0, 1; t25 * t49 + t72, -t25 * t47 + t49 * t78, t24, t25 * pkin(2) + t24 * qJ(3) + t70; t23 * t49 - t47 * t77, -t23 * t47 - t49 * t77, t22, t23 * pkin(2) + t22 * qJ(3) + t65; t32 * t49 + t79, -t32 * t47 + t50 * t49, -t30, t32 * pkin(2) - t30 * qJ(3) + t71; 0, 0, 0, 1; t12, -t11, t24, t63; t10, -t9, t22, t59; t17, -t16, -t30, t64; 0, 0, 0, 1; t4, -t3, t11, t61; t2, -t1, t9, t58; t6, -t5, t16, t62; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(5) + t3 * qJ(6) + t61; t2, t9, t1, t2 * pkin(5) + t1 * qJ(6) + t58; t6, t16, t5, t6 * pkin(5) + t5 * qJ(6) + t62; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
