% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:05
% EndTime: 2019-12-05 17:23:05
% DurationCPUTime: 0.24s
% Computational Cost: add. (315->63), mult. (825->93), div. (0->0), fcn. (1108->14), ass. (0->53)
t41 = sin(pkin(6));
t44 = cos(pkin(6));
t45 = cos(pkin(5));
t42 = sin(pkin(5));
t52 = cos(qJ(2));
t72 = t42 * t52;
t25 = -t41 * t72 + t45 * t44;
t40 = sin(pkin(11));
t43 = cos(pkin(11));
t49 = sin(qJ(2));
t68 = t45 * t52;
t28 = -t40 * t68 - t43 * t49;
t74 = t40 * t42;
t19 = -t28 * t41 + t44 * t74;
t76 = cos(qJ(3));
t73 = t42 * t49;
t71 = t43 * t42;
t69 = t45 * t49;
t65 = qJ(1) + 0;
t64 = t41 * t76;
t63 = t44 * t76;
t62 = t43 * pkin(1) + pkin(7) * t74 + 0;
t61 = t45 * pkin(7) + t65;
t60 = t42 * t64;
t26 = -t40 * t49 + t43 * t68;
t18 = -t26 * t41 - t44 * t71;
t59 = t40 * pkin(1) - pkin(7) * t71 + 0;
t29 = -t40 * t69 + t43 * t52;
t58 = t29 * pkin(2) + t19 * pkin(8) + t62;
t57 = pkin(2) * t73 + t25 * pkin(8) + t61;
t48 = sin(qJ(3));
t10 = t29 * t76 + (t28 * t44 + t41 * t74) * t48;
t9 = -t28 * t63 + t29 * t48 - t40 * t60;
t56 = t10 * pkin(3) + t9 * pkin(9) + t58;
t16 = -t45 * t64 + t48 * t73 - t63 * t72;
t17 = t45 * t41 * t48 + (t44 * t48 * t52 + t76 * t49) * t42;
t55 = t17 * pkin(3) + t16 * pkin(9) + t57;
t27 = t40 * t52 + t43 * t69;
t54 = t27 * pkin(2) + t18 * pkin(8) + t59;
t7 = -t26 * t63 + t27 * t48 + t43 * t60;
t8 = t27 * t76 + (t26 * t44 - t41 * t71) * t48;
t53 = t8 * pkin(3) + t7 * pkin(9) + t54;
t51 = cos(qJ(4));
t50 = cos(qJ(5));
t47 = sin(qJ(4));
t46 = sin(qJ(5));
t12 = t17 * t51 + t25 * t47;
t11 = t17 * t47 - t25 * t51;
t4 = t10 * t51 + t19 * t47;
t3 = t10 * t47 - t19 * t51;
t2 = t18 * t47 + t8 * t51;
t1 = -t18 * t51 + t8 * t47;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t43, -t40, 0, 0; t40, t43, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t29, t28, t74, t62; t27, t26, -t71, t59; t73, t72, t45, t61; 0, 0, 0, 1; t10, -t9, t19, t58; t8, -t7, t18, t54; t17, -t16, t25, t57; 0, 0, 0, 1; t4, -t3, t9, t56; t2, -t1, t7, t53; t12, -t11, t16, t55; 0, 0, 0, 1; t4 * t50 + t9 * t46, -t4 * t46 + t9 * t50, t3, t4 * pkin(4) + t3 * pkin(10) + t56; t2 * t50 + t7 * t46, -t2 * t46 + t7 * t50, t1, t2 * pkin(4) + t1 * pkin(10) + t53; t12 * t50 + t16 * t46, -t12 * t46 + t16 * t50, t11, t12 * pkin(4) + t11 * pkin(10) + t55; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
