% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:23:22
% EndTime: 2018-11-23 17:23:22
% DurationCPUTime: 0.23s
% Computational Cost: add. (717->88), mult. (534->103), div. (0->0), fcn. (579->22), ass. (0->64)
t46 = pkin(6) - qJ(2);
t31 = cos(t46) / 0.2e1;
t45 = pkin(6) + qJ(2);
t40 = cos(t45);
t81 = t31 - t40 / 0.2e1;
t30 = sin(t45) / 0.2e1;
t38 = sin(t46);
t20 = t30 - t38 / 0.2e1;
t48 = sin(pkin(6));
t54 = sin(qJ(1));
t32 = t54 * t48;
t58 = cos(qJ(1));
t78 = t58 * t48;
t77 = pkin(7) + 0;
t44 = qJ(2) + pkin(12);
t51 = sin(qJ(5));
t76 = pkin(5) * t51 + pkin(9);
t50 = pkin(8) + qJ(3);
t14 = pkin(2) * t20 - t48 * t50;
t57 = cos(qJ(2));
t35 = pkin(2) * t57 + pkin(1);
t75 = t14 * t58 + t35 * t54 + 0;
t74 = pkin(6) - t44;
t73 = pkin(6) + t44;
t64 = sin(t73) / 0.2e1;
t69 = sin(t74);
t17 = t64 - t69 / 0.2e1;
t39 = cos(t44);
t8 = t17 * t58 + t39 * t54;
t72 = pkin(3) * t8 + t75;
t71 = -t54 * t14 + t35 * t58 + 0;
t70 = cos(t73);
t49 = cos(pkin(6));
t68 = pkin(2) * t81 + t49 * t50 + t77;
t10 = -t17 * t54 + t39 * t58;
t67 = pkin(3) * t10 + t71;
t65 = cos(t74) / 0.2e1;
t19 = t65 - t70 / 0.2e1;
t66 = pkin(3) * t19 + t68;
t36 = sin(t44);
t60 = t70 / 0.2e1 + t65;
t7 = t36 * t54 - t58 * t60;
t63 = pkin(9) * t7 + t72;
t9 = t36 * t58 + t54 * t60;
t62 = pkin(9) * t9 + t67;
t18 = t69 / 0.2e1 + t64;
t61 = -pkin(9) * t18 + t66;
t59 = -pkin(11) - pkin(10);
t56 = cos(qJ(4));
t55 = cos(qJ(5));
t53 = sin(qJ(2));
t52 = sin(qJ(4));
t47 = qJ(5) + qJ(6);
t43 = cos(t47);
t42 = sin(t47);
t34 = pkin(5) * t55 + pkin(4);
t21 = t31 + t40 / 0.2e1;
t12 = t19 * t56 + t49 * t52;
t11 = t19 * t52 - t49 * t56;
t4 = t10 * t56 + t32 * t52;
t3 = t10 * t52 - t32 * t56;
t2 = -t52 * t78 + t56 * t8;
t1 = t52 * t8 + t56 * t78;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t58, -t54, 0, 0; t54, t58, 0, 0; 0, 0, 1, t77; 0, 0, 0, 1; -t20 * t54 + t57 * t58, -t21 * t54 - t53 * t58, t32, pkin(1) * t58 + pkin(8) * t32 + 0; t20 * t58 + t54 * t57, t21 * t58 - t53 * t54, -t78, pkin(1) * t54 - pkin(8) * t78 + 0; t81, t30 + t38 / 0.2e1, t49, pkin(8) * t49 + t77; 0, 0, 0, 1; t10, -t9, t32, t71; t8, -t7, -t78, t75; t19, t18, t49, t68; 0, 0, 0, 1; t4, -t3, t9, t62; t2, -t1, t7, t63; t12, -t11, -t18, t61; 0, 0, 0, 1; t4 * t55 + t51 * t9, -t4 * t51 + t55 * t9, t3, pkin(4) * t4 + pkin(10) * t3 + t62; t2 * t55 + t51 * t7, -t2 * t51 + t55 * t7, t1, pkin(4) * t2 + pkin(10) * t1 + t63; t12 * t55 - t18 * t51, -t12 * t51 - t18 * t55, t11, pkin(4) * t12 + pkin(10) * t11 + t61; 0, 0, 0, 1; t4 * t43 + t42 * t9, -t4 * t42 + t43 * t9, t3, -t3 * t59 + t4 * t34 + t76 * t9 + t67; t2 * t43 + t42 * t7, -t2 * t42 + t43 * t7, t1, -t1 * t59 + t2 * t34 + t7 * t76 + t72; t12 * t43 - t18 * t42, -t12 * t42 - t18 * t43, t11, -t11 * t59 + t12 * t34 - t18 * t76 + t66; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
