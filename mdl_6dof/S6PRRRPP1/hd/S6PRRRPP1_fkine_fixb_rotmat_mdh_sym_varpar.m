% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2018-11-23 15:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRRPP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:20:01
% EndTime: 2018-11-23 15:20:02
% DurationCPUTime: 0.19s
% Computational Cost: add. (656->75), mult. (703->85), div. (0->0), fcn. (793->16), ass. (0->57)
t48 = sin(pkin(10));
t50 = cos(pkin(10));
t55 = sin(qJ(2));
t75 = pkin(6) - qJ(2);
t68 = cos(t75) / 0.2e1;
t74 = pkin(6) + qJ(2);
t70 = cos(t74);
t60 = t68 + t70 / 0.2e1;
t23 = t48 * t55 - t50 * t60;
t53 = sin(qJ(4));
t80 = t23 * t53;
t25 = t48 * t60 + t50 * t55;
t79 = t25 * t53;
t67 = sin(t74) / 0.2e1;
t69 = sin(t75);
t33 = t67 + t69 / 0.2e1;
t78 = t33 * t53;
t49 = sin(pkin(6));
t77 = t48 * t49;
t76 = t50 * t49;
t73 = qJ(1) + 0;
t72 = t50 * pkin(1) + pkin(7) * t77 + 0;
t51 = cos(pkin(6));
t71 = t51 * pkin(7) + t73;
t66 = t48 * pkin(1) - pkin(7) * t76 + 0;
t34 = t67 - t69 / 0.2e1;
t58 = cos(qJ(2));
t26 = -t48 * t34 + t50 * t58;
t65 = t26 * pkin(2) + t25 * pkin(8) + t72;
t35 = t68 - t70 / 0.2e1;
t64 = t35 * pkin(2) - t33 * pkin(8) + t71;
t24 = t50 * t34 + t48 * t58;
t63 = t24 * pkin(2) + t23 * pkin(8) + t66;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t13 = t26 * t54 - t57 * t77;
t14 = t26 * t57 + t54 * t77;
t56 = cos(qJ(4));
t41 = t56 * pkin(4) + pkin(3);
t52 = -qJ(5) - pkin(9);
t62 = pkin(4) * t79 - t13 * t52 + t14 * t41 + t65;
t27 = t35 * t54 - t51 * t57;
t28 = t35 * t57 + t51 * t54;
t61 = -pkin(4) * t78 - t27 * t52 + t28 * t41 + t64;
t11 = t24 * t54 + t57 * t76;
t12 = t24 * t57 - t54 * t76;
t59 = pkin(4) * t80 - t11 * t52 + t12 * t41 + t63;
t47 = qJ(4) + pkin(11);
t43 = cos(t47);
t42 = sin(t47);
t6 = t28 * t43 - t33 * t42;
t5 = t28 * t42 + t33 * t43;
t4 = t14 * t43 + t25 * t42;
t3 = t14 * t42 - t25 * t43;
t2 = t12 * t43 + t23 * t42;
t1 = t12 * t42 - t23 * t43;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t50, -t48, 0, 0; t48, t50, 0, 0; 0, 0, 1, t73; 0, 0, 0, 1; t26, -t25, t77, t72; t24, -t23, -t76, t66; t35, t33, t51, t71; 0, 0, 0, 1; t14, -t13, t25, t65; t12, -t11, t23, t63; t28, -t27, -t33, t64; 0, 0, 0, 1; t14 * t56 + t79, -t14 * t53 + t25 * t56, t13, t14 * pkin(3) + t13 * pkin(9) + t65; t12 * t56 + t80, -t12 * t53 + t23 * t56, t11, t12 * pkin(3) + t11 * pkin(9) + t63; t28 * t56 - t78, -t28 * t53 - t33 * t56, t27, t28 * pkin(3) + t27 * pkin(9) + t64; 0, 0, 0, 1; t4, -t3, t13, t62; t2, -t1, t11, t59; t6, -t5, t27, t61; 0, 0, 0, 1; t4, t13, t3, t4 * pkin(5) + t3 * qJ(6) + t62; t2, t11, t1, t2 * pkin(5) + t1 * qJ(6) + t59; t6, t27, t5, t6 * pkin(5) + t5 * qJ(6) + t61; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
