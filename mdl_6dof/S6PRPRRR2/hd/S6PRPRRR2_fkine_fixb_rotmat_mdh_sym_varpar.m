% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2018-11-23 15:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRPRRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:03:42
% EndTime: 2018-11-23 15:03:42
% DurationCPUTime: 0.25s
% Computational Cost: add. (717->88), mult. (534->103), div. (0->0), fcn. (579->22), ass. (0->64)
t46 = pkin(6) - qJ(2);
t32 = cos(t46) / 0.2e1;
t45 = pkin(6) + qJ(2);
t40 = cos(t45);
t81 = t32 - t40 / 0.2e1;
t31 = sin(t45) / 0.2e1;
t38 = sin(t46);
t20 = t31 - t38 / 0.2e1;
t48 = sin(pkin(11));
t49 = sin(pkin(6));
t29 = t48 * t49;
t50 = cos(pkin(11));
t78 = t50 * t49;
t44 = qJ(2) + pkin(12);
t77 = qJ(1) + 0;
t53 = sin(qJ(5));
t76 = pkin(5) * t53 + pkin(8);
t52 = pkin(7) + qJ(3);
t14 = t20 * pkin(2) - t49 * t52;
t58 = cos(qJ(2));
t35 = t58 * pkin(2) + pkin(1);
t75 = t50 * t14 + t48 * t35 + 0;
t74 = pkin(6) - t44;
t73 = pkin(6) + t44;
t65 = sin(t73) / 0.2e1;
t69 = sin(t74);
t17 = t65 - t69 / 0.2e1;
t39 = cos(t44);
t8 = t50 * t17 + t48 * t39;
t72 = t8 * pkin(3) + t75;
t71 = -t48 * t14 + t50 * t35 + 0;
t70 = cos(t73);
t51 = cos(pkin(6));
t68 = t81 * pkin(2) + t51 * t52 + t77;
t10 = -t48 * t17 + t50 * t39;
t67 = t10 * pkin(3) + t71;
t66 = cos(t74) / 0.2e1;
t19 = t66 - t70 / 0.2e1;
t64 = t19 * pkin(3) + t68;
t36 = sin(t44);
t60 = t70 / 0.2e1 + t66;
t7 = t48 * t36 - t50 * t60;
t63 = t7 * pkin(8) + t72;
t9 = t50 * t36 + t48 * t60;
t62 = t9 * pkin(8) + t67;
t18 = t69 / 0.2e1 + t65;
t61 = -t18 * pkin(8) + t64;
t59 = -pkin(10) - pkin(9);
t57 = cos(qJ(4));
t56 = cos(qJ(5));
t55 = sin(qJ(2));
t54 = sin(qJ(4));
t47 = qJ(5) + qJ(6);
t43 = cos(t47);
t42 = sin(t47);
t34 = t56 * pkin(5) + pkin(4);
t21 = t32 + t40 / 0.2e1;
t12 = t19 * t57 + t51 * t54;
t11 = t19 * t54 - t51 * t57;
t4 = t10 * t57 + t54 * t29;
t3 = t10 * t54 - t57 * t29;
t2 = -t54 * t78 + t8 * t57;
t1 = t8 * t54 + t57 * t78;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t50, -t48, 0, 0; t48, t50, 0, 0; 0, 0, 1, t77; 0, 0, 0, 1; -t48 * t20 + t50 * t58, -t48 * t21 - t50 * t55, t29, t50 * pkin(1) + pkin(7) * t29 + 0; t50 * t20 + t48 * t58, t50 * t21 - t48 * t55, -t78, t48 * pkin(1) - pkin(7) * t78 + 0; t81, t31 + t38 / 0.2e1, t51, t51 * pkin(7) + t77; 0, 0, 0, 1; t10, -t9, t29, t71; t8, -t7, -t78, t75; t19, t18, t51, t68; 0, 0, 0, 1; t4, -t3, t9, t62; t2, -t1, t7, t63; t12, -t11, -t18, t61; 0, 0, 0, 1; t4 * t56 + t9 * t53, -t4 * t53 + t9 * t56, t3, t4 * pkin(4) + t3 * pkin(9) + t62; t2 * t56 + t7 * t53, -t2 * t53 + t7 * t56, t1, t2 * pkin(4) + t1 * pkin(9) + t63; t12 * t56 - t18 * t53, -t12 * t53 - t18 * t56, t11, t12 * pkin(4) + t11 * pkin(9) + t61; 0, 0, 0, 1; t4 * t43 + t9 * t42, -t4 * t42 + t9 * t43, t3, -t3 * t59 + t4 * t34 + t76 * t9 + t67; t2 * t43 + t7 * t42, -t2 * t42 + t7 * t43, t1, -t1 * t59 + t2 * t34 + t76 * t7 + t72; t12 * t43 - t18 * t42, -t12 * t42 - t18 * t43, t11, -t11 * t59 + t12 * t34 - t76 * t18 + t64; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
