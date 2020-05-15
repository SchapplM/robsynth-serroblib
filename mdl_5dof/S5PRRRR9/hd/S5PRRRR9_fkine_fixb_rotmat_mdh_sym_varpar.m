% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:18:53
% EndTime: 2019-12-05 17:18:53
% DurationCPUTime: 0.15s
% Computational Cost: add. (179->61), mult. (371->79), div. (0->0), fcn. (511->12), ass. (0->44)
t25 = sin(pkin(10));
t26 = sin(pkin(5));
t53 = t25 * t26;
t31 = sin(qJ(2));
t52 = t26 * t31;
t33 = cos(qJ(3));
t51 = t26 * t33;
t34 = cos(qJ(2));
t50 = t26 * t34;
t27 = cos(pkin(10));
t49 = t27 * t26;
t28 = cos(pkin(5));
t48 = t28 * t31;
t47 = t28 * t34;
t46 = qJ(1) + 0;
t29 = sin(qJ(4));
t45 = pkin(4) * t29 + pkin(7);
t44 = t27 * pkin(1) + pkin(6) * t53 + 0;
t43 = t28 * pkin(6) + t46;
t10 = -t25 * t48 + t27 * t34;
t42 = t10 * pkin(2) + t44;
t41 = pkin(2) * t52 + t43;
t40 = t25 * pkin(1) - pkin(6) * t49 + 0;
t9 = t25 * t47 + t27 * t31;
t39 = t9 * pkin(7) + t42;
t8 = t25 * t34 + t27 * t48;
t38 = t8 * pkin(2) + t40;
t37 = -pkin(7) * t50 + t41;
t7 = t25 * t31 - t27 * t47;
t36 = t7 * pkin(7) + t38;
t35 = -pkin(9) - pkin(8);
t32 = cos(qJ(4));
t30 = sin(qJ(3));
t24 = qJ(4) + qJ(5);
t20 = cos(t24);
t19 = sin(t24);
t18 = t32 * pkin(4) + pkin(3);
t12 = t28 * t30 + t31 * t51;
t11 = -t28 * t33 + t30 * t52;
t4 = t10 * t33 + t30 * t53;
t3 = t10 * t30 - t25 * t51;
t2 = -t30 * t49 + t8 * t33;
t1 = t8 * t30 + t33 * t49;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t25, 0, 0; t25, t27, 0, 0; 0, 0, 1, t46; 0, 0, 0, 1; t10, -t9, t53, t44; t8, -t7, -t49, t40; t52, t50, t28, t43; 0, 0, 0, 1; t4, -t3, t9, t39; t2, -t1, t7, t36; t12, -t11, -t50, t37; 0, 0, 0, 1; t9 * t29 + t4 * t32, -t4 * t29 + t9 * t32, t3, t4 * pkin(3) + t3 * pkin(8) + t39; t2 * t32 + t7 * t29, -t2 * t29 + t7 * t32, t1, t2 * pkin(3) + t1 * pkin(8) + t36; t12 * t32 - t29 * t50, -t12 * t29 - t32 * t50, t11, t12 * pkin(3) + t11 * pkin(8) + t37; 0, 0, 0, 1; t9 * t19 + t4 * t20, -t4 * t19 + t9 * t20, t3, t4 * t18 - t3 * t35 + t45 * t9 + t42; t7 * t19 + t2 * t20, -t2 * t19 + t7 * t20, t1, -t1 * t35 + t2 * t18 + t45 * t7 + t38; t12 * t20 - t19 * t50, -t12 * t19 - t20 * t50, t11, -t11 * t35 + t12 * t18 - t45 * t50 + t41; 0, 0, 0, 1;];
T_ges = t5;
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
