% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2018-11-23 15:13
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRPRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:13:06
% EndTime: 2018-11-23 15:13:06
% DurationCPUTime: 0.20s
% Computational Cost: add. (572->74), mult. (646->75), div. (0->0), fcn. (727->14), ass. (0->60)
t76 = pkin(4) + pkin(8);
t37 = sin(pkin(10));
t39 = cos(pkin(10));
t43 = sin(qJ(2));
t67 = pkin(6) - qJ(2);
t56 = cos(t67) / 0.2e1;
t66 = pkin(6) + qJ(2);
t60 = cos(t66);
t47 = t56 + t60 / 0.2e1;
t16 = t37 * t43 - t39 * t47;
t75 = t16 * pkin(8);
t18 = t37 * t47 + t39 * t43;
t74 = t18 * pkin(8);
t55 = sin(t66) / 0.2e1;
t59 = sin(t67);
t24 = t55 + t59 / 0.2e1;
t73 = t24 * pkin(8);
t44 = cos(qJ(5));
t72 = t44 * pkin(5) + t76;
t71 = cos(qJ(3));
t38 = sin(pkin(6));
t70 = t37 * t38;
t69 = t39 * t38;
t68 = cos(pkin(6));
t65 = qJ(1) + 0;
t64 = t38 * t71;
t41 = sin(qJ(5));
t63 = pkin(5) * t41 + qJ(4);
t62 = t39 * pkin(1) + pkin(7) * t70 + 0;
t61 = t68 * pkin(7) + t65;
t25 = t55 - t59 / 0.2e1;
t45 = cos(qJ(2));
t19 = -t37 * t25 + t39 * t45;
t58 = t19 * pkin(2) + t62;
t26 = t56 - t60 / 0.2e1;
t57 = t26 * pkin(2) + t61;
t42 = sin(qJ(3));
t12 = t19 * t71 + t42 * t70;
t54 = t12 * pkin(3) + t58;
t21 = t26 * t71 + t68 * t42;
t53 = t21 * pkin(3) + t57;
t52 = t37 * pkin(1) - pkin(7) * t69 + 0;
t17 = t39 * t25 + t37 * t45;
t51 = t17 * pkin(2) + t52;
t10 = t17 * t71 - t42 * t69;
t50 = t10 * pkin(3) + t51;
t11 = t19 * t42 - t37 * t64;
t49 = t11 * qJ(4) + t54;
t20 = t26 * t42 - t68 * t71;
t48 = t20 * qJ(4) + t53;
t9 = t17 * t42 + t39 * t64;
t46 = t9 * qJ(4) + t50;
t40 = -qJ(6) - pkin(9);
t6 = t20 * t41 - t24 * t44;
t5 = t20 * t44 + t24 * t41;
t4 = t11 * t41 + t18 * t44;
t3 = t11 * t44 - t18 * t41;
t2 = t16 * t44 + t9 * t41;
t1 = -t16 * t41 + t9 * t44;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t39, -t37, 0, 0; t37, t39, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t19, -t18, t70, t62; t17, -t16, -t69, t52; t26, t24, t68, t61; 0, 0, 0, 1; t12, -t11, t18, t58 + t74; t10, -t9, t16, t51 + t75; t21, -t20, -t24, t57 - t73; 0, 0, 0, 1; t18, -t12, t11, t49 + t74; t16, -t10, t9, t46 + t75; -t24, -t21, t20, t48 - t73; 0, 0, 0, 1; t4, t3, t12, t12 * pkin(9) + t76 * t18 + t49; t2, t1, t10, t10 * pkin(9) + t76 * t16 + t46; t6, t5, t21, t21 * pkin(9) - t76 * t24 + t48; 0, 0, 0, 1; t4, t3, t12, t63 * t11 - t12 * t40 + t72 * t18 + t54; t2, t1, t10, -t10 * t40 + t72 * t16 + t63 * t9 + t50; t6, t5, t21, t63 * t20 - t21 * t40 - t72 * t24 + t53; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
