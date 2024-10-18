% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RRRRR14_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:10
% EndTime: 2024-09-27 18:42:10
% DurationCPUTime: 0.02s
% Computational Cost: add. (208->61), mult. (119->71), div. (0->0), fcn. (154->22), ass. (0->48)
t43 = pkin(8) + pkin(9);
t41 = cos(qJ(3));
t22 = t41 * pkin(3) + pkin(2);
t35 = qJ(1) + qJ(2);
t26 = sin(t35);
t36 = sin(pkin(5));
t10 = t26 * t36;
t28 = cos(t35);
t52 = t28 * t36;
t39 = sin(qJ(3));
t51 = t36 * t39;
t37 = cos(pkin(5));
t50 = t37 * t39;
t49 = t37 * t41;
t34 = qJ(3) + qJ(4);
t48 = pkin(6) + 0;
t40 = sin(qJ(1));
t47 = t40 * pkin(1) + 0;
t42 = cos(qJ(1));
t46 = t42 * pkin(1) + 0;
t45 = pkin(7) + t48;
t24 = pkin(5) - t34;
t23 = pkin(5) + t34;
t44 = pkin(4) * sin(qJ(4)) * t41 + (cos(qJ(4)) * pkin(4) + pkin(3)) * t39;
t33 = pkin(10) + t43;
t29 = qJ(5) + t34;
t27 = cos(t34);
t25 = sin(t34);
t20 = cos(t29);
t19 = sin(t29);
t18 = -qJ(5) + t24;
t17 = qJ(5) + t23;
t16 = cos(t23);
t15 = sin(t24);
t14 = cos(t17);
t13 = sin(t18);
t12 = cos(t24) / 0.2e1;
t11 = sin(t23) / 0.2e1;
t9 = cos(t18) / 0.2e1;
t8 = sin(t17) / 0.2e1;
t7 = pkin(4) * t27 + t22;
t6 = pkin(3) * t50 - t36 * t43;
t5 = t12 + t16 / 0.2e1;
t4 = t11 - t15 / 0.2e1;
t3 = t9 + t14 / 0.2e1;
t2 = t8 - t13 / 0.2e1;
t1 = -t36 * t33 + t44 * t37;
t21 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t42, -t40, 0, 0; t40, t42, 0, 0; 0, 0, 1, t48; t28, -t26, 0, t46; t26, t28, 0, t47; 0, 0, 1, t45; -t26 * t50 + t28 * t41, -t26 * t49 - t28 * t39, t10, t28 * pkin(2) + pkin(8) * t10 + t46; t26 * t41 + t28 * t50, -t26 * t39 + t28 * t49, -t52, t26 * pkin(2) - pkin(8) * t52 + t47; t51, t36 * t41, t37, t37 * pkin(8) + t45; -t26 * t4 + t28 * t27, -t28 * t25 - t26 * t5, t10, t28 * t22 - t26 * t6 + t46; t26 * t27 + t28 * t4, -t26 * t25 + t28 * t5, -t52, t26 * t22 + t28 * t6 + t47; t12 - t16 / 0.2e1, t11 + t15 / 0.2e1, t37, pkin(3) * t51 + t37 * t43 + t45; -t26 * t2 + t28 * t20, -t28 * t19 - t26 * t3, t10, -t26 * t1 + t28 * t7 + t46; t28 * t2 + t26 * t20, -t26 * t19 + t28 * t3, -t52, t28 * t1 + t26 * t7 + t47; t9 - t14 / 0.2e1, t8 + t13 / 0.2e1, t37, t37 * t33 + t44 * t36 + t45;];
Tc_stack = t21;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
