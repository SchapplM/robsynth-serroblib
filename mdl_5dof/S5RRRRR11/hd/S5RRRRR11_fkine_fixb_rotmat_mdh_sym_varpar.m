% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRRR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:39
% EndTime: 2019-12-31 22:38:39
% DurationCPUTime: 0.14s
% Computational Cost: add. (179->61), mult. (371->76), div. (0->0), fcn. (511->12), ass. (0->45)
t25 = sin(pkin(5));
t29 = sin(qJ(2));
t54 = t25 * t29;
t33 = cos(qJ(2));
t53 = t25 * t33;
t30 = sin(qJ(1));
t52 = t30 * t25;
t51 = t30 * t29;
t50 = t30 * t33;
t34 = cos(qJ(1));
t49 = t34 * t25;
t48 = t34 * t29;
t47 = t34 * t33;
t46 = pkin(6) + 0;
t27 = sin(qJ(4));
t45 = pkin(4) * t27 + pkin(8);
t26 = cos(pkin(5));
t44 = t26 * pkin(7) + t46;
t43 = t34 * pkin(1) + pkin(7) * t52 + 0;
t42 = pkin(2) * t54 + t44;
t12 = -t26 * t51 + t47;
t41 = t12 * pkin(2) + t43;
t40 = t30 * pkin(1) - pkin(7) * t49 + 0;
t10 = t26 * t48 + t50;
t39 = t10 * pkin(2) + t40;
t11 = t26 * t50 + t48;
t38 = t11 * pkin(8) + t41;
t37 = -pkin(8) * t53 + t42;
t9 = -t26 * t47 + t51;
t36 = t9 * pkin(8) + t39;
t35 = -pkin(10) - pkin(9);
t32 = cos(qJ(3));
t31 = cos(qJ(4));
t28 = sin(qJ(3));
t24 = qJ(4) + qJ(5);
t20 = cos(t24);
t19 = sin(t24);
t18 = t31 * pkin(4) + pkin(3);
t8 = t26 * t28 + t32 * t54;
t7 = -t26 * t32 + t28 * t54;
t4 = t12 * t32 + t28 * t52;
t3 = t12 * t28 - t32 * t52;
t2 = t10 * t32 - t28 * t49;
t1 = t10 * t28 + t32 * t49;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t30, 0, 0; t30, t34, 0, 0; 0, 0, 1, t46; 0, 0, 0, 1; t12, -t11, t52, t43; t10, -t9, -t49, t40; t54, t53, t26, t44; 0, 0, 0, 1; t4, -t3, t11, t38; t2, -t1, t9, t36; t8, -t7, -t53, t37; 0, 0, 0, 1; t11 * t27 + t4 * t31, t11 * t31 - t4 * t27, t3, t4 * pkin(3) + t3 * pkin(9) + t38; t2 * t31 + t9 * t27, -t2 * t27 + t9 * t31, t1, t2 * pkin(3) + t1 * pkin(9) + t36; -t27 * t53 + t8 * t31, -t8 * t27 - t31 * t53, t7, t8 * pkin(3) + t7 * pkin(9) + t37; 0, 0, 0, 1; t11 * t19 + t4 * t20, t11 * t20 - t4 * t19, t3, t45 * t11 + t4 * t18 - t3 * t35 + t41; t9 * t19 + t2 * t20, -t2 * t19 + t9 * t20, t1, -t1 * t35 + t2 * t18 + t45 * t9 + t39; -t19 * t53 + t8 * t20, -t8 * t19 - t20 * t53, t7, t8 * t18 - t7 * t35 - t45 * t53 + t42; 0, 0, 0, 1;];
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
