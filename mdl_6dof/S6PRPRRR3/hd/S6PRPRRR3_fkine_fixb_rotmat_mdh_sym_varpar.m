% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR3
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

function T_c_mdh = S6PRPRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:04:26
% EndTime: 2018-11-23 15:04:26
% DurationCPUTime: 0.21s
% Computational Cost: add. (541->89), mult. (472->103), div. (0->0), fcn. (534->18), ass. (0->55)
t42 = sin(pkin(12));
t72 = t42 * pkin(3);
t45 = cos(pkin(12));
t32 = t45 * pkin(3) + pkin(2);
t43 = sin(pkin(11));
t44 = sin(pkin(6));
t71 = t43 * t44;
t46 = cos(pkin(11));
t70 = t46 * t44;
t47 = cos(pkin(6));
t69 = t47 * t42;
t48 = -pkin(8) - qJ(3);
t68 = pkin(6) - qJ(2);
t67 = pkin(6) + qJ(2);
t41 = pkin(12) + qJ(4);
t66 = t43 * pkin(1) + 0;
t65 = t42 * t71;
t64 = qJ(1) + 0;
t63 = t46 * pkin(1) + pkin(7) * t71 + 0;
t62 = t47 * pkin(7) + t64;
t61 = cos(t67);
t60 = sin(t68);
t59 = cos(t68) / 0.2e1;
t58 = sin(t67) / 0.2e1;
t57 = -pkin(7) * t70 + t66;
t50 = sin(qJ(2));
t53 = t59 + t61 / 0.2e1;
t14 = t43 * t53 + t46 * t50;
t23 = t58 - t60 / 0.2e1;
t52 = cos(qJ(2));
t15 = -t43 * t23 + t46 * t52;
t34 = cos(t41);
t21 = pkin(4) * t34 + t32;
t33 = sin(t41);
t25 = pkin(4) * t33 + t72;
t40 = -pkin(9) + t48;
t56 = -t14 * t40 + t15 * t21 + t25 * t71 + t63;
t22 = t58 + t60 / 0.2e1;
t24 = t59 - t61 / 0.2e1;
t55 = t24 * t21 + t22 * t40 + t47 * t25 + t62;
t12 = t43 * t50 - t46 * t53;
t13 = t46 * t23 + t43 * t52;
t54 = t13 * t21 - t12 * t40 + (-pkin(7) - t25) * t70 + t66;
t51 = cos(qJ(6));
t49 = sin(qJ(6));
t35 = qJ(5) + t41;
t31 = cos(t35);
t30 = sin(t35);
t8 = t24 * t31 + t47 * t30;
t7 = t24 * t30 - t47 * t31;
t4 = t15 * t31 + t30 * t71;
t3 = t15 * t30 - t31 * t71;
t2 = t13 * t31 - t30 * t70;
t1 = t13 * t30 + t31 * t70;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t46, -t43, 0, 0; t43, t46, 0, 0; 0, 0, 1, t64; 0, 0, 0, 1; t15, -t14, t71, t63; t13, -t12, -t70, t57; t24, t22, t47, t62; 0, 0, 0, 1; t15 * t45 + t65, -t15 * t42 + t45 * t71, t14, t15 * pkin(2) + t14 * qJ(3) + t63; t13 * t45 - t42 * t70, -t13 * t42 - t45 * t70, t12, t13 * pkin(2) + t12 * qJ(3) + t57; t24 * t45 + t69, -t24 * t42 + t47 * t45, -t22, t24 * pkin(2) - t22 * qJ(3) + t62; 0, 0, 0, 1; t15 * t34 + t33 * t71, -t15 * t33 + t34 * t71, t14, pkin(3) * t65 - t14 * t48 + t15 * t32 + t63; t13 * t34 - t33 * t70, -t13 * t33 - t34 * t70, t12, -t12 * t48 + t13 * t32 + (-pkin(7) - t72) * t70 + t66; t24 * t34 + t47 * t33, -t24 * t33 + t47 * t34, -t22, pkin(3) * t69 + t22 * t48 + t24 * t32 + t62; 0, 0, 0, 1; t4, -t3, t14, t56; t2, -t1, t12, t54; t8, -t7, -t22, t55; 0, 0, 0, 1; t14 * t49 + t4 * t51, t14 * t51 - t4 * t49, t3, t4 * pkin(5) + t3 * pkin(10) + t56; t12 * t49 + t2 * t51, t12 * t51 - t2 * t49, t1, t2 * pkin(5) + t1 * pkin(10) + t54; -t22 * t49 + t8 * t51, -t22 * t51 - t8 * t49, t7, t8 * pkin(5) + t7 * pkin(10) + t55; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
