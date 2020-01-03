% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRR14_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:43
% EndTime: 2019-12-31 19:16:43
% DurationCPUTime: 0.22s
% Computational Cost: add. (315->63), mult. (825->92), div. (0->0), fcn. (1108->14), ass. (0->56)
t45 = cos(pkin(5));
t40 = sin(pkin(11));
t52 = cos(qJ(1));
t71 = t52 * t40;
t43 = cos(pkin(11));
t49 = sin(qJ(1));
t72 = t49 * t43;
t28 = -t45 * t72 - t71;
t41 = sin(pkin(6));
t44 = cos(pkin(6));
t42 = sin(pkin(5));
t73 = t49 * t42;
t19 = -t28 * t41 + t44 * t73;
t76 = t42 * t43;
t25 = -t41 * t76 + t45 * t44;
t79 = cos(qJ(3));
t77 = t42 * t40;
t74 = t49 * t40;
t70 = t52 * t42;
t69 = t52 * t43;
t68 = qJ(2) * t42;
t67 = pkin(7) + 0;
t64 = t41 * t79;
t63 = t44 * t79;
t62 = t45 * qJ(2) + t67;
t61 = t52 * pkin(1) + t49 * t68 + 0;
t60 = t42 * t64;
t26 = t45 * t69 - t74;
t18 = -t26 * t41 - t44 * t70;
t59 = t49 * pkin(1) - t52 * t68 + 0;
t29 = -t45 * t74 + t69;
t58 = t29 * pkin(2) + t19 * pkin(8) + t61;
t57 = pkin(2) * t77 + t25 * pkin(8) + t62;
t48 = sin(qJ(3));
t11 = -t28 * t63 + t29 * t48 - t49 * t60;
t12 = t29 * t79 + (t28 * t44 + t41 * t73) * t48;
t56 = t12 * pkin(3) + t11 * pkin(9) + t58;
t16 = -t45 * t64 + t48 * t77 - t63 * t76;
t17 = t45 * t41 * t48 + (t43 * t44 * t48 + t79 * t40) * t42;
t55 = t17 * pkin(3) + t16 * pkin(9) + t57;
t27 = t45 * t71 + t72;
t54 = t27 * pkin(2) + t18 * pkin(8) + t59;
t10 = t27 * t79 + (t26 * t44 - t41 * t70) * t48;
t9 = -t26 * t63 + t27 * t48 + t52 * t60;
t53 = t10 * pkin(3) + t9 * pkin(9) + t54;
t51 = cos(qJ(4));
t50 = cos(qJ(5));
t47 = sin(qJ(4));
t46 = sin(qJ(5));
t8 = t17 * t51 + t25 * t47;
t7 = t17 * t47 - t25 * t51;
t4 = t12 * t51 + t19 * t47;
t3 = t12 * t47 - t19 * t51;
t2 = t10 * t51 + t18 * t47;
t1 = t10 * t47 - t18 * t51;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t52, -t49, 0, 0; t49, t52, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t29, t28, t73, t61; t27, t26, -t70, t59; t77, t76, t45, t62; 0, 0, 0, 1; t12, -t11, t19, t58; t10, -t9, t18, t54; t17, -t16, t25, t57; 0, 0, 0, 1; t4, -t3, t11, t56; t2, -t1, t9, t53; t8, -t7, t16, t55; 0, 0, 0, 1; t11 * t46 + t4 * t50, t11 * t50 - t4 * t46, t3, t4 * pkin(4) + t3 * pkin(10) + t56; t2 * t50 + t9 * t46, -t2 * t46 + t9 * t50, t1, t2 * pkin(4) + t1 * pkin(10) + t53; t16 * t46 + t8 * t50, t16 * t50 - t8 * t46, t7, t8 * pkin(4) + t7 * pkin(10) + t55; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
