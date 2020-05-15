% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRPR13_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:09
% EndTime: 2019-12-31 21:43:10
% DurationCPUTime: 0.15s
% Computational Cost: add. (157->56), mult. (349->59), div. (0->0), fcn. (481->10), ass. (0->43)
t56 = pkin(4) + pkin(8);
t26 = cos(pkin(5));
t33 = cos(qJ(2));
t34 = cos(qJ(1));
t46 = t34 * t33;
t29 = sin(qJ(2));
t30 = sin(qJ(1));
t50 = t30 * t29;
t12 = -t26 * t46 + t50;
t55 = t12 * pkin(8);
t47 = t34 * t29;
t49 = t30 * t33;
t14 = t26 * t49 + t47;
t54 = t14 * pkin(8);
t25 = sin(pkin(5));
t53 = t25 * t29;
t52 = t25 * t33;
t51 = t30 * t25;
t48 = t34 * t25;
t45 = pkin(6) + 0;
t44 = pkin(8) * t52;
t43 = t26 * pkin(7) + t45;
t42 = t34 * pkin(1) + pkin(7) * t51 + 0;
t41 = pkin(2) * t53 + t43;
t15 = -t26 * t50 + t46;
t40 = t15 * pkin(2) + t42;
t39 = t30 * pkin(1) - pkin(7) * t48 + 0;
t13 = t26 * t47 + t49;
t38 = t13 * pkin(2) + t39;
t28 = sin(qJ(3));
t32 = cos(qJ(3));
t5 = t15 * t28 - t32 * t51;
t6 = t15 * t32 + t28 * t51;
t37 = t6 * pkin(3) + t5 * qJ(4) + t40;
t10 = -t26 * t32 + t28 * t53;
t11 = t26 * t28 + t32 * t53;
t36 = t11 * pkin(3) + t10 * qJ(4) + t41;
t3 = t13 * t28 + t32 * t48;
t4 = t13 * t32 - t28 * t48;
t35 = t4 * pkin(3) + t3 * qJ(4) + t38;
t31 = cos(qJ(5));
t27 = sin(qJ(5));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t30, 0, 0; t30, t34, 0, 0; 0, 0, 1, t45; 0, 0, 0, 1; t15, -t14, t51, t42; t13, -t12, -t48, t39; t53, t52, t26, t43; 0, 0, 0, 1; t6, -t5, t14, t40 + t54; t4, -t3, t12, t38 + t55; t11, -t10, -t52, t41 - t44; 0, 0, 0, 1; t14, -t6, t5, t37 + t54; t12, -t4, t3, t35 + t55; -t52, -t11, t10, t36 - t44; 0, 0, 0, 1; t14 * t31 + t5 * t27, -t14 * t27 + t5 * t31, t6, t6 * pkin(9) + t14 * t56 + t37; t12 * t31 + t3 * t27, -t12 * t27 + t3 * t31, t4, t4 * pkin(9) + t12 * t56 + t35; t10 * t27 - t31 * t52, t10 * t31 + t27 * t52, t11, t11 * pkin(9) - t52 * t56 + t36; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
