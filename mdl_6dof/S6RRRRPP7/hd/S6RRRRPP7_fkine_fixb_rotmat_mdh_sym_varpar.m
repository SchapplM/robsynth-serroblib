% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2018-11-23 18:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRPP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:08:46
% EndTime: 2018-11-23 18:08:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (656->75), mult. (703->85), div. (0->0), fcn. (793->16), ass. (0->57)
t53 = sin(qJ(2));
t54 = sin(qJ(1));
t58 = cos(qJ(1));
t75 = pkin(6) - qJ(2);
t68 = cos(t75) / 0.2e1;
t74 = pkin(6) + qJ(2);
t70 = cos(t74);
t60 = t68 + t70 / 0.2e1;
t25 = t54 * t53 - t58 * t60;
t51 = sin(qJ(4));
t80 = t25 * t51;
t27 = t58 * t53 + t54 * t60;
t79 = t27 * t51;
t67 = sin(t74) / 0.2e1;
t69 = sin(t75);
t33 = t67 + t69 / 0.2e1;
t78 = t33 * t51;
t48 = sin(pkin(6));
t77 = t54 * t48;
t76 = t58 * t48;
t73 = pkin(7) + 0;
t49 = cos(pkin(6));
t72 = t49 * pkin(8) + t73;
t71 = t58 * pkin(1) + pkin(8) * t77 + 0;
t66 = t54 * pkin(1) - pkin(8) * t76 + 0;
t35 = t68 - t70 / 0.2e1;
t65 = t35 * pkin(2) - t33 * pkin(9) + t72;
t34 = t67 - t69 / 0.2e1;
t57 = cos(qJ(2));
t28 = -t54 * t34 + t58 * t57;
t64 = t28 * pkin(2) + t27 * pkin(9) + t71;
t26 = t58 * t34 + t54 * t57;
t63 = t26 * pkin(2) + t25 * pkin(9) + t66;
t52 = sin(qJ(3));
t56 = cos(qJ(3));
t23 = t35 * t52 - t49 * t56;
t24 = t35 * t56 + t49 * t52;
t55 = cos(qJ(4));
t41 = t55 * pkin(4) + pkin(3);
t50 = -qJ(5) - pkin(10);
t62 = -pkin(4) * t78 - t23 * t50 + t24 * t41 + t65;
t13 = t28 * t52 - t56 * t77;
t14 = t28 * t56 + t52 * t77;
t61 = pkin(4) * t79 - t13 * t50 + t14 * t41 + t64;
t11 = t26 * t52 + t56 * t76;
t12 = t26 * t56 - t52 * t76;
t59 = pkin(4) * t80 - t11 * t50 + t12 * t41 + t63;
t47 = qJ(4) + pkin(11);
t43 = cos(t47);
t42 = sin(t47);
t6 = t24 * t43 - t33 * t42;
t5 = t24 * t42 + t33 * t43;
t4 = t14 * t43 + t27 * t42;
t3 = t14 * t42 - t27 * t43;
t2 = t12 * t43 + t25 * t42;
t1 = t12 * t42 - t25 * t43;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t58, -t54, 0, 0; t54, t58, 0, 0; 0, 0, 1, t73; 0, 0, 0, 1; t28, -t27, t77, t71; t26, -t25, -t76, t66; t35, t33, t49, t72; 0, 0, 0, 1; t14, -t13, t27, t64; t12, -t11, t25, t63; t24, -t23, -t33, t65; 0, 0, 0, 1; t14 * t55 + t79, -t14 * t51 + t27 * t55, t13, t14 * pkin(3) + t13 * pkin(10) + t64; t12 * t55 + t80, -t12 * t51 + t25 * t55, t11, t12 * pkin(3) + t11 * pkin(10) + t63; t24 * t55 - t78, -t24 * t51 - t33 * t55, t23, t24 * pkin(3) + t23 * pkin(10) + t65; 0, 0, 0, 1; t4, -t3, t13, t61; t2, -t1, t11, t59; t6, -t5, t23, t62; 0, 0, 0, 1; t4, t13, t3, t4 * pkin(5) + t3 * qJ(6) + t61; t2, t11, t1, t2 * pkin(5) + t1 * qJ(6) + t59; t6, t23, t5, t6 * pkin(5) + t5 * qJ(6) + t62; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
