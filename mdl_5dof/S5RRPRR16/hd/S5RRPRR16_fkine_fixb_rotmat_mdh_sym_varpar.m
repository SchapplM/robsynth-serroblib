% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-29 19:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR16_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 19:29:14
% EndTime: 2019-12-29 19:29:15
% DurationCPUTime: 0.31s
% Computational Cost: add. (137->54), mult. (294->59), div. (0->0), fcn. (405->10), ass. (0->40)
t26 = sin(pkin(5));
t30 = sin(qJ(2));
t21 = t26 * t30;
t34 = cos(qJ(2));
t52 = t26 * t34;
t31 = sin(qJ(1));
t20 = t31 * t26;
t51 = t31 * t30;
t50 = t31 * t34;
t35 = cos(qJ(1));
t49 = t35 * t26;
t48 = t35 * t30;
t47 = t35 * t34;
t46 = pkin(6) + 0;
t45 = pkin(7) * t49;
t44 = t31 * pkin(1) + 0;
t27 = cos(pkin(5));
t43 = t27 * pkin(7) + t46;
t42 = t35 * pkin(1) + pkin(7) * t20 + 0;
t10 = -t27 * t47 + t51;
t11 = t27 * t48 + t50;
t41 = t11 * pkin(2) + t10 * qJ(3) + t44;
t12 = t27 * t50 + t48;
t13 = -t27 * t51 + t47;
t40 = t13 * pkin(2) + t12 * qJ(3) + t42;
t39 = pkin(2) * t21 - qJ(3) * t52 + t43;
t38 = t27 * pkin(3) + pkin(8) * t21 + t39;
t37 = pkin(3) * t20 + t13 * pkin(8) + t40;
t36 = t11 * pkin(8) + (-pkin(3) - pkin(7)) * t49 + t41;
t33 = cos(qJ(4));
t32 = cos(qJ(5));
t29 = sin(qJ(4));
t28 = sin(qJ(5));
t9 = t27 * t33 - t29 * t52;
t8 = t27 * t29 + t33 * t52;
t4 = t10 * t29 - t33 * t49;
t3 = t10 * t33 + t29 * t49;
t2 = t12 * t29 + t33 * t20;
t1 = -t12 * t33 + t29 * t20;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t31, 0, 0; t31, t35, 0, 0; 0, 0, 1, t46; 0, 0, 0, 1; t13, -t12, t20, t42; t11, -t10, -t49, t44 - t45; t21, t52, t27, t43; 0, 0, 0, 1; t20, -t13, t12, t40; -t49, -t11, t10, t41 - t45; t27, -t21, -t52, t39; 0, 0, 0, 1; t2, -t1, t13, t37; t4, t3, t11, t36; t9, -t8, t21, t38; 0, 0, 0, 1; t13 * t28 + t2 * t32, t13 * t32 - t2 * t28, t1, t2 * pkin(4) + t1 * pkin(9) + t37; t11 * t28 + t4 * t32, t11 * t32 - t4 * t28, -t3, t4 * pkin(4) - t3 * pkin(9) + t36; t28 * t21 + t9 * t32, t32 * t21 - t9 * t28, t8, t9 * pkin(4) + t8 * pkin(9) + t38; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
