% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2018-11-23 14:57
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:57:46
% EndTime: 2018-11-23 14:57:46
% DurationCPUTime: 0.16s
% Computational Cost: add. (552->78), mult. (536->86), div. (0->0), fcn. (603->16), ass. (0->51)
t40 = sin(pkin(10));
t41 = sin(pkin(6));
t71 = t40 * t41;
t43 = cos(pkin(10));
t70 = t43 * t41;
t39 = sin(pkin(11));
t44 = cos(pkin(6));
t69 = t44 * t39;
t68 = pkin(6) - qJ(2);
t67 = pkin(6) + qJ(2);
t66 = t40 * pkin(1) + 0;
t65 = t39 * t71;
t64 = qJ(1) + 0;
t63 = t43 * pkin(1) + pkin(7) * t71 + 0;
t62 = t44 * pkin(7) + t64;
t61 = cos(t67);
t60 = sin(t68);
t59 = cos(t68) / 0.2e1;
t58 = sin(t67) / 0.2e1;
t57 = -pkin(7) * t70 + t66;
t47 = sin(qJ(2));
t53 = t59 + t61 / 0.2e1;
t16 = t40 * t53 + t43 * t47;
t22 = t58 - t60 / 0.2e1;
t49 = cos(qJ(2));
t17 = -t40 * t22 + t43 * t49;
t42 = cos(pkin(11));
t32 = t42 * pkin(3) + pkin(2);
t45 = -pkin(8) - qJ(3);
t56 = pkin(3) * t65 - t16 * t45 + t17 * t32 + t63;
t21 = t58 + t60 / 0.2e1;
t23 = t59 - t61 / 0.2e1;
t55 = pkin(3) * t69 + t21 * t45 + t23 * t32 + t62;
t38 = pkin(11) + qJ(4);
t33 = sin(t38);
t34 = cos(t38);
t5 = t17 * t33 - t34 * t71;
t6 = t17 * t34 + t33 * t71;
t54 = t6 * pkin(4) + t5 * qJ(5) + t56;
t10 = t23 * t33 - t44 * t34;
t11 = t23 * t34 + t44 * t33;
t52 = t11 * pkin(4) + t10 * qJ(5) + t55;
t14 = t40 * t47 - t43 * t53;
t15 = t43 * t22 + t40 * t49;
t51 = -t14 * t45 + t15 * t32 + (-pkin(3) * t39 - pkin(7)) * t70 + t66;
t3 = t15 * t33 + t34 * t70;
t4 = t15 * t34 - t33 * t70;
t50 = t4 * pkin(4) + t3 * qJ(5) + t51;
t48 = cos(qJ(6));
t46 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t43, -t40, 0, 0; t40, t43, 0, 0; 0, 0, 1, t64; 0, 0, 0, 1; t17, -t16, t71, t63; t15, -t14, -t70, t57; t23, t21, t44, t62; 0, 0, 0, 1; t17 * t42 + t65, -t17 * t39 + t42 * t71, t16, t17 * pkin(2) + t16 * qJ(3) + t63; t15 * t42 - t39 * t70, -t15 * t39 - t42 * t70, t14, t15 * pkin(2) + t14 * qJ(3) + t57; t23 * t42 + t69, -t23 * t39 + t44 * t42, -t21, t23 * pkin(2) - t21 * qJ(3) + t62; 0, 0, 0, 1; t6, -t5, t16, t56; t4, -t3, t14, t51; t11, -t10, -t21, t55; 0, 0, 0, 1; t16, -t6, t5, t54; t14, -t4, t3, t50; -t21, -t11, t10, t52; 0, 0, 0, 1; t16 * t48 + t5 * t46, -t16 * t46 + t5 * t48, t6, t16 * pkin(5) + t6 * pkin(9) + t54; t14 * t48 + t3 * t46, -t14 * t46 + t3 * t48, t4, t14 * pkin(5) + t4 * pkin(9) + t50; t10 * t46 - t21 * t48, t10 * t48 + t21 * t46, t11, -t21 * pkin(5) + t11 * pkin(9) + t52; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
