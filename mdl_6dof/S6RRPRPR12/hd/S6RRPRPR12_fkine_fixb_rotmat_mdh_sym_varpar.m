% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRPR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:08:17
% EndTime: 2018-11-23 17:08:17
% DurationCPUTime: 0.18s
% Computational Cost: add. (505->79), mult. (495->85), div. (0->0), fcn. (547->16), ass. (0->52)
t45 = sin(qJ(2));
t46 = sin(qJ(1));
t50 = cos(qJ(1));
t69 = pkin(6) - qJ(2);
t60 = cos(t69) / 0.2e1;
t68 = pkin(6) + qJ(2);
t62 = cos(t68);
t52 = t60 + t62 / 0.2e1;
t14 = t46 * t45 - t50 * t52;
t44 = sin(qJ(4));
t73 = t14 * t44;
t16 = t50 * t45 + t46 * t52;
t72 = t16 * t44;
t59 = sin(t68) / 0.2e1;
t61 = sin(t69);
t25 = t59 + t61 / 0.2e1;
t71 = t25 * t44;
t40 = sin(pkin(6));
t32 = t46 * t40;
t70 = t50 * t40;
t67 = pkin(7) + 0;
t66 = pkin(8) * t70;
t65 = t46 * pkin(1) + 0;
t41 = cos(pkin(6));
t64 = t41 * pkin(8) + t67;
t63 = t50 * pkin(1) + pkin(8) * t32 + 0;
t58 = t59 - t61 / 0.2e1;
t49 = cos(qJ(2));
t15 = t46 * t49 + t50 * t58;
t57 = t15 * pkin(2) + t14 * qJ(3) + t65;
t26 = t60 - t62 / 0.2e1;
t56 = t26 * pkin(2) - t25 * qJ(3) + t64;
t17 = -t46 * t58 + t50 * t49;
t55 = t17 * pkin(2) + t16 * qJ(3) + t63;
t48 = cos(qJ(4));
t33 = t48 * pkin(4) + pkin(3);
t42 = -qJ(5) - pkin(9);
t54 = -pkin(4) * t71 - t26 * t42 + t41 * t33 + t56;
t53 = pkin(4) * t72 - t17 * t42 + t33 * t32 + t55;
t51 = -t15 * t42 + pkin(4) * t73 + (-pkin(8) - t33) * t70 + t57;
t47 = cos(qJ(6));
t43 = sin(qJ(6));
t39 = qJ(4) + pkin(11);
t35 = cos(t39);
t34 = sin(t39);
t7 = -t25 * t34 + t41 * t35;
t6 = t25 * t35 + t41 * t34;
t4 = t14 * t34 - t35 * t70;
t3 = t14 * t35 + t34 * t70;
t2 = t16 * t34 + t35 * t32;
t1 = -t16 * t35 + t34 * t32;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t50, -t46, 0, 0; t46, t50, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t17, -t16, t32, t63; t15, -t14, -t70, t65 - t66; t26, t25, t41, t64; 0, 0, 0, 1; t32, -t17, t16, t55; -t70, -t15, t14, t57 - t66; t41, -t26, -t25, t56; 0, 0, 0, 1; t48 * t32 + t72, t16 * t48 - t44 * t32, t17, pkin(3) * t32 + t17 * pkin(9) + t55; -t48 * t70 + t73, t14 * t48 + t44 * t70, t15, t15 * pkin(9) + (-pkin(3) - pkin(8)) * t70 + t57; t41 * t48 - t71, -t25 * t48 - t41 * t44, t26, t41 * pkin(3) + t26 * pkin(9) + t56; 0, 0, 0, 1; t2, -t1, t17, t53; t4, t3, t15, t51; t7, -t6, t26, t54; 0, 0, 0, 1; t17 * t43 + t2 * t47, t17 * t47 - t2 * t43, t1, t2 * pkin(5) + t1 * pkin(10) + t53; t15 * t43 + t4 * t47, t15 * t47 - t4 * t43, -t3, t4 * pkin(5) - t3 * pkin(10) + t51; t26 * t43 + t7 * t47, t26 * t47 - t7 * t43, t6, t7 * pkin(5) + t6 * pkin(10) + t54; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
