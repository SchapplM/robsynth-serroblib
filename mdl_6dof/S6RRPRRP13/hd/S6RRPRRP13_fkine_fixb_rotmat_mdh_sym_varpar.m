% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRP13_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:19:31
% EndTime: 2018-11-23 17:19:31
% DurationCPUTime: 0.18s
% Computational Cost: add. (508->73), mult. (558->74), div. (0->0), fcn. (622->14), ass. (0->56)
t36 = sin(pkin(6));
t42 = sin(qJ(1));
t30 = t42 * t36;
t46 = cos(qJ(1));
t70 = t46 * t36;
t69 = pkin(6) - qJ(2);
t68 = pkin(6) + qJ(2);
t67 = pkin(7) + 0;
t66 = pkin(8) * t70;
t65 = t42 * pkin(1) + 0;
t39 = sin(qJ(5));
t64 = pkin(5) * t39 + pkin(9);
t37 = cos(pkin(6));
t63 = t37 * pkin(8) + t67;
t62 = t46 * pkin(1) + pkin(8) * t30 + 0;
t61 = cos(t68);
t60 = sin(t69);
t59 = cos(t69) / 0.2e1;
t58 = sin(t68) / 0.2e1;
t57 = t59 + t61 / 0.2e1;
t41 = sin(qJ(2));
t16 = t42 * t41 - t46 * t57;
t45 = cos(qJ(2));
t51 = t58 - t60 / 0.2e1;
t17 = t42 * t45 + t46 * t51;
t56 = t17 * pkin(2) + t16 * qJ(3) + t65;
t24 = t58 + t60 / 0.2e1;
t25 = t59 - t61 / 0.2e1;
t55 = t25 * pkin(2) - t24 * qJ(3) + t63;
t18 = t46 * t41 + t42 * t57;
t19 = -t42 * t51 + t46 * t45;
t54 = t19 * pkin(2) + t18 * qJ(3) + t62;
t53 = t37 * pkin(3) + t55;
t52 = pkin(3) * t30 + t54;
t50 = t25 * pkin(9) + t53;
t49 = t19 * pkin(9) + t52;
t48 = (-pkin(3) - pkin(8)) * t70 + t56;
t47 = t17 * pkin(9) + t48;
t44 = cos(qJ(4));
t43 = cos(qJ(5));
t40 = sin(qJ(4));
t38 = -qJ(6) - pkin(10);
t31 = t43 * pkin(5) + pkin(4);
t15 = -t24 * t40 + t37 * t44;
t14 = t24 * t44 + t37 * t40;
t10 = t16 * t40 - t44 * t70;
t9 = t16 * t44 + t40 * t70;
t8 = t18 * t40 + t44 * t30;
t7 = -t18 * t44 + t40 * t30;
t6 = t15 * t43 + t25 * t39;
t5 = -t15 * t39 + t25 * t43;
t4 = t10 * t43 + t17 * t39;
t3 = -t10 * t39 + t17 * t43;
t2 = t19 * t39 + t8 * t43;
t1 = t19 * t43 - t8 * t39;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t46, -t42, 0, 0; t42, t46, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t19, -t18, t30, t62; t17, -t16, -t70, t65 - t66; t25, t24, t37, t63; 0, 0, 0, 1; t30, -t19, t18, t54; -t70, -t17, t16, t56 - t66; t37, -t25, -t24, t55; 0, 0, 0, 1; t8, -t7, t19, t49; t10, t9, t17, t47; t15, -t14, t25, t50; 0, 0, 0, 1; t2, t1, t7, t8 * pkin(4) + t7 * pkin(10) + t49; t4, t3, -t9, t10 * pkin(4) - t9 * pkin(10) + t47; t6, t5, t14, t15 * pkin(4) + t14 * pkin(10) + t50; 0, 0, 0, 1; t2, t1, t7, t64 * t19 + t8 * t31 - t7 * t38 + t52; t4, t3, -t9, t10 * t31 + t64 * t17 + t9 * t38 + t48; t6, t5, t14, -t14 * t38 + t15 * t31 + t64 * t25 + t53; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
