% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRPP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:11:51
% EndTime: 2019-12-05 16:11:52
% DurationCPUTime: 0.12s
% Computational Cost: add. (109->44), mult. (216->50), div. (0->0), fcn. (301->8), ass. (0->36)
t29 = sin(pkin(7));
t33 = sin(qJ(2));
t50 = t29 * t33;
t35 = cos(qJ(2));
t49 = t29 * t35;
t31 = cos(pkin(7));
t48 = t31 * t33;
t47 = t31 * t35;
t32 = sin(qJ(3));
t46 = t32 * t35;
t22 = t33 * t32;
t34 = cos(qJ(3));
t45 = t33 * t34;
t44 = t34 * t35;
t27 = qJ(1) + 0;
t43 = t31 * pkin(1) + t29 * pkin(5) + 0;
t42 = t29 * pkin(1) - t31 * pkin(5) + 0;
t41 = pkin(2) * t47 + pkin(6) * t48 + t43;
t40 = t33 * pkin(2) - t35 * pkin(6) + t27;
t39 = pkin(2) * t49 + pkin(6) * t50 + t42;
t38 = pkin(3) * t45 + qJ(4) * t22 + t40;
t11 = -t29 * t34 + t31 * t46;
t12 = t29 * t32 + t31 * t44;
t37 = t12 * pkin(3) + t11 * qJ(4) + t41;
t7 = t29 * t46 + t31 * t34;
t8 = t29 * t44 - t31 * t32;
t36 = t8 * pkin(3) + t7 * qJ(4) + t39;
t30 = cos(pkin(8));
t28 = sin(pkin(8));
t10 = -t35 * t28 + t30 * t45;
t9 = t28 * t45 + t35 * t30;
t4 = t12 * t30 + t28 * t48;
t3 = t12 * t28 - t30 * t48;
t2 = t28 * t50 + t8 * t30;
t1 = t8 * t28 - t30 * t50;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t31, -t29, 0, 0; t29, t31, 0, 0; 0, 0, 1, t27; 0, 0, 0, 1; t47, -t48, t29, t43; t49, -t50, -t31, t42; t33, t35, 0, t27; 0, 0, 0, 1; t12, -t11, t48, t41; t8, -t7, t50, t39; t45, -t22, -t35, t40; 0, 0, 0, 1; t4, -t3, t11, t37; t2, -t1, t7, t36; t10, -t9, t22, t38; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(4) + t3 * qJ(5) + t37; t2, t7, t1, t2 * pkin(4) + t1 * qJ(5) + t36; t10, t22, t9, t10 * pkin(4) + t9 * qJ(5) + t38; 0, 0, 0, 1;];
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
