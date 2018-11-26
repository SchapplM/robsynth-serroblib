% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRP3
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

function T_c_mdh = S6PRRPRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:12:13
% EndTime: 2018-11-23 15:12:13
% DurationCPUTime: 0.19s
% Computational Cost: add. (656->75), mult. (703->85), div. (0->0), fcn. (793->16), ass. (0->57)
t49 = sin(pkin(10));
t52 = cos(pkin(10));
t56 = sin(qJ(2));
t75 = pkin(6) - qJ(2);
t68 = cos(t75) / 0.2e1;
t74 = pkin(6) + qJ(2);
t70 = cos(t74);
t60 = t68 + t70 / 0.2e1;
t23 = t49 * t56 - t52 * t60;
t48 = sin(pkin(11));
t80 = t23 * t48;
t25 = t49 * t60 + t52 * t56;
t79 = t25 * t48;
t67 = sin(t74) / 0.2e1;
t69 = sin(t75);
t33 = t67 + t69 / 0.2e1;
t78 = t33 * t48;
t50 = sin(pkin(6));
t77 = t49 * t50;
t76 = t52 * t50;
t73 = qJ(1) + 0;
t72 = t52 * pkin(1) + pkin(7) * t77 + 0;
t53 = cos(pkin(6));
t71 = t53 * pkin(7) + t73;
t66 = t49 * pkin(1) - pkin(7) * t76 + 0;
t34 = t67 - t69 / 0.2e1;
t58 = cos(qJ(2));
t26 = -t49 * t34 + t52 * t58;
t65 = t26 * pkin(2) + t25 * pkin(8) + t72;
t35 = t68 - t70 / 0.2e1;
t64 = t35 * pkin(2) - t33 * pkin(8) + t71;
t24 = t52 * t34 + t49 * t58;
t63 = t24 * pkin(2) + t23 * pkin(8) + t66;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t13 = t26 * t55 - t57 * t77;
t14 = t26 * t57 + t55 * t77;
t51 = cos(pkin(11));
t41 = t51 * pkin(4) + pkin(3);
t54 = -pkin(9) - qJ(4);
t62 = pkin(4) * t79 - t13 * t54 + t14 * t41 + t65;
t27 = t35 * t55 - t53 * t57;
t28 = t35 * t57 + t53 * t55;
t61 = -pkin(4) * t78 - t27 * t54 + t28 * t41 + t64;
t11 = t24 * t55 + t57 * t76;
t12 = t24 * t57 - t55 * t76;
t59 = pkin(4) * t80 - t11 * t54 + t12 * t41 + t63;
t47 = pkin(11) + qJ(5);
t43 = cos(t47);
t42 = sin(t47);
t6 = t28 * t43 - t33 * t42;
t5 = t28 * t42 + t33 * t43;
t4 = t14 * t43 + t25 * t42;
t3 = t14 * t42 - t25 * t43;
t2 = t12 * t43 + t23 * t42;
t1 = t12 * t42 - t23 * t43;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t52, -t49, 0, 0; t49, t52, 0, 0; 0, 0, 1, t73; 0, 0, 0, 1; t26, -t25, t77, t72; t24, -t23, -t76, t66; t35, t33, t53, t71; 0, 0, 0, 1; t14, -t13, t25, t65; t12, -t11, t23, t63; t28, -t27, -t33, t64; 0, 0, 0, 1; t14 * t51 + t79, -t14 * t48 + t25 * t51, t13, t14 * pkin(3) + t13 * qJ(4) + t65; t12 * t51 + t80, -t12 * t48 + t23 * t51, t11, t12 * pkin(3) + t11 * qJ(4) + t63; t28 * t51 - t78, -t28 * t48 - t33 * t51, t27, t28 * pkin(3) + t27 * qJ(4) + t64; 0, 0, 0, 1; t4, -t3, t13, t62; t2, -t1, t11, t59; t6, -t5, t27, t61; 0, 0, 0, 1; t4, t13, t3, t4 * pkin(5) + t3 * qJ(6) + t62; t2, t11, t1, t2 * pkin(5) + t1 * qJ(6) + t59; t6, t27, t5, t6 * pkin(5) + t5 * qJ(6) + t61; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
