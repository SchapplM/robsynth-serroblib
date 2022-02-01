% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)
% T_c_stack [(5+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:58:56
% EndTime: 2022-01-23 08:58:56
% DurationCPUTime: 0.18s
% Computational Cost: add. (122->68), mult. (200->86), div. (0->0), fcn. (280->10), ass. (0->48)
t25 = sin(pkin(7));
t50 = t25 * qJ(3) + pkin(1);
t24 = sin(pkin(8));
t49 = t24 * qJ(4) + pkin(2);
t29 = sin(qJ(5));
t48 = t24 * t29;
t31 = cos(qJ(5));
t47 = t24 * t31;
t46 = t25 * t24;
t27 = cos(pkin(8));
t45 = t25 * t27;
t26 = cos(pkin(9));
t28 = cos(pkin(7));
t44 = t28 * t26;
t30 = sin(qJ(1));
t43 = t30 * t25;
t42 = t30 * t28;
t32 = cos(qJ(1));
t41 = t32 * t24;
t40 = t32 * t25;
t39 = t32 * t28;
t22 = pkin(5) + 0;
t38 = t30 * qJ(2) + 0;
t37 = qJ(4) * t27 - qJ(2);
t36 = -t32 * qJ(2) + 0;
t35 = -t28 * qJ(3) + t22;
t23 = sin(pkin(9));
t34 = t23 * pkin(4) - t26 * pkin(6) + qJ(3);
t18 = t26 * pkin(4) + t23 * pkin(6) + pkin(3);
t33 = -t18 * t24 + t37;
t17 = t27 * pkin(3) + t49;
t16 = pkin(2) * t28 + t50;
t15 = -t24 * pkin(3) + t37;
t14 = t30 * t24 + t27 * t39;
t13 = t24 * t39 - t30 * t27;
t12 = t27 * t42 - t41;
t11 = t24 * t42 + t32 * t27;
t10 = t26 * t48 + t31 * t27;
t9 = t25 * t23 + t27 * t44;
t8 = -t28 * t23 + t26 * t45;
t7 = t23 * t45 + t44;
t6 = t18 * t27 + t49;
t5 = t17 * t28 + t50;
t4 = t14 * t23 - t26 * t40;
t3 = t12 * t23 - t26 * t43;
t2 = t28 * t47 - t9 * t29;
t1 = t34 * t25 + t6 * t28 + pkin(1);
t19 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t32, -t30, 0, 0; t30, t32, 0, 0; 0, 0, 1, t22; t39, -t40, t30, t32 * pkin(1) + t38; t42, -t43, -t32, t30 * pkin(1) + t36; t25, t28, 0, t22; t14, -t13, t40, t16 * t32 + t38; t12, -t11, t43, t16 * t30 + t36; t45, -t46, -t28, t25 * pkin(2) + t35; t14 * t26 + t23 * t40, -t4, t13, -t15 * t30 + t5 * t32 + 0; t12 * t26 + t23 * t43, -t3, t11, t15 * t32 + t5 * t30 + 0; t8, -t7, t46, t17 * t25 + t35; (t28 * t48 + t9 * t31) * t32 + t30 * (t26 * t47 - t29 * t27), -t30 * t10 + t2 * t32, t4, t1 * t32 - t33 * t30 + 0; (-t26 * t41 + t9 * t30) * t31 + t11 * t29, t32 * t10 + t2 * t30, t3, t1 * t30 + t33 * t32 + 0; t29 * t46 + t8 * t31, -t8 * t29 + t31 * t46, t7, t6 * t25 - t34 * t28 + t22;];
Tc_stack = t19;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
