% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:08:28
% EndTime: 2018-11-23 15:08:28
% DurationCPUTime: 0.17s
% Computational Cost: add. (552->78), mult. (536->86), div. (0->0), fcn. (603->16), ass. (0->51)
t39 = sin(pkin(10));
t40 = sin(pkin(6));
t71 = t39 * t40;
t41 = cos(pkin(10));
t70 = t41 * t40;
t42 = cos(pkin(6));
t45 = sin(qJ(3));
t69 = t42 * t45;
t68 = pkin(6) - qJ(2);
t67 = pkin(6) + qJ(2);
t66 = t39 * pkin(1) + 0;
t65 = t45 * t71;
t64 = qJ(1) + 0;
t63 = t41 * pkin(1) + pkin(7) * t71 + 0;
t62 = t42 * pkin(7) + t64;
t61 = cos(t67);
t60 = sin(t68);
t59 = cos(t68) / 0.2e1;
t58 = sin(t67) / 0.2e1;
t57 = -pkin(7) * t70 + t66;
t46 = sin(qJ(2));
t53 = t59 + t61 / 0.2e1;
t16 = t39 * t53 + t41 * t46;
t22 = t58 - t60 / 0.2e1;
t49 = cos(qJ(2));
t17 = -t39 * t22 + t41 * t49;
t48 = cos(qJ(3));
t32 = t48 * pkin(3) + pkin(2);
t43 = -qJ(4) - pkin(8);
t56 = pkin(3) * t65 - t16 * t43 + t17 * t32 + t63;
t21 = t58 + t60 / 0.2e1;
t23 = t59 - t61 / 0.2e1;
t55 = pkin(3) * t69 + t21 * t43 + t23 * t32 + t62;
t38 = qJ(3) + pkin(11);
t33 = sin(t38);
t34 = cos(t38);
t5 = t17 * t33 - t34 * t71;
t6 = t17 * t34 + t33 * t71;
t54 = t6 * pkin(4) + t5 * qJ(5) + t56;
t10 = t23 * t33 - t42 * t34;
t11 = t23 * t34 + t42 * t33;
t52 = t11 * pkin(4) + t10 * qJ(5) + t55;
t14 = t39 * t46 - t41 * t53;
t15 = t41 * t22 + t39 * t49;
t51 = -t14 * t43 + t15 * t32 + (-pkin(3) * t45 - pkin(7)) * t70 + t66;
t3 = t15 * t33 + t34 * t70;
t4 = t15 * t34 - t33 * t70;
t50 = t4 * pkin(4) + t3 * qJ(5) + t51;
t47 = cos(qJ(6));
t44 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t41, -t39, 0, 0; t39, t41, 0, 0; 0, 0, 1, t64; 0, 0, 0, 1; t17, -t16, t71, t63; t15, -t14, -t70, t57; t23, t21, t42, t62; 0, 0, 0, 1; t17 * t48 + t65, -t17 * t45 + t48 * t71, t16, t17 * pkin(2) + t16 * pkin(8) + t63; t15 * t48 - t45 * t70, -t15 * t45 - t48 * t70, t14, t15 * pkin(2) + t14 * pkin(8) + t57; t23 * t48 + t69, -t23 * t45 + t42 * t48, -t21, t23 * pkin(2) - t21 * pkin(8) + t62; 0, 0, 0, 1; t6, -t5, t16, t56; t4, -t3, t14, t51; t11, -t10, -t21, t55; 0, 0, 0, 1; t16, -t6, t5, t54; t14, -t4, t3, t50; -t21, -t11, t10, t52; 0, 0, 0, 1; t16 * t47 + t5 * t44, -t16 * t44 + t5 * t47, t6, t16 * pkin(5) + t6 * pkin(9) + t54; t14 * t47 + t3 * t44, -t14 * t44 + t3 * t47, t4, t14 * pkin(5) + t4 * pkin(9) + t50; t10 * t44 - t21 * t47, t10 * t47 + t21 * t44, t11, -t21 * pkin(5) + t11 * pkin(9) + t52; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
