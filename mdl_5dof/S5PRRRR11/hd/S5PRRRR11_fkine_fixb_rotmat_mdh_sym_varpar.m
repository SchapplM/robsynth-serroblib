% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5PRRRR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:44:59
% EndTime: 2024-09-27 21:44:59
% DurationCPUTime: 0.00s
% Computational Cost: add. (208->61), mult. (119->71), div. (0->0), fcn. (154->22), ass. (0->48)
t43 = pkin(7) + pkin(8);
t42 = cos(qJ(3));
t22 = t42 * pkin(3) + pkin(2);
t33 = pkin(10) + qJ(2);
t23 = sin(t33);
t37 = sin(pkin(5));
t10 = t23 * t37;
t24 = cos(t33);
t52 = t24 * t37;
t41 = sin(qJ(3));
t51 = t37 * t41;
t39 = cos(pkin(5));
t50 = t39 * t41;
t49 = t39 * t42;
t35 = qJ(3) + qJ(4);
t36 = sin(pkin(10));
t48 = t36 * pkin(1) + 0;
t38 = cos(pkin(10));
t47 = t38 * pkin(1) + 0;
t46 = qJ(1) + 0;
t26 = pkin(5) - t35;
t25 = pkin(5) + t35;
t45 = pkin(6) + t46;
t44 = pkin(4) * sin(qJ(4)) * t42 + (cos(qJ(4)) * pkin(4) + pkin(3)) * t41;
t34 = pkin(9) + t43;
t31 = qJ(5) + t35;
t28 = cos(t35);
t27 = sin(t35);
t20 = cos(t31);
t19 = sin(t31);
t18 = -qJ(5) + t26;
t17 = qJ(5) + t25;
t16 = cos(t25);
t15 = sin(t26);
t14 = cos(t17);
t13 = sin(t18);
t12 = cos(t26) / 0.2e1;
t11 = sin(t25) / 0.2e1;
t9 = cos(t18) / 0.2e1;
t8 = sin(t17) / 0.2e1;
t7 = pkin(4) * t28 + t22;
t6 = pkin(3) * t50 - t37 * t43;
t5 = t12 + t16 / 0.2e1;
t4 = t11 - t15 / 0.2e1;
t3 = t9 + t14 / 0.2e1;
t2 = t8 - t13 / 0.2e1;
t1 = -t37 * t34 + t44 * t39;
t21 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t38, -t36, 0, 0; t36, t38, 0, 0; 0, 0, 1, t46; t24, -t23, 0, t47; t23, t24, 0, t48; 0, 0, 1, t45; -t23 * t50 + t24 * t42, -t23 * t49 - t24 * t41, t10, t24 * pkin(2) + pkin(7) * t10 + t47; t23 * t42 + t24 * t50, -t23 * t41 + t24 * t49, -t52, t23 * pkin(2) - pkin(7) * t52 + t48; t51, t37 * t42, t39, t39 * pkin(7) + t45; -t23 * t4 + t24 * t28, -t23 * t5 - t24 * t27, t10, t24 * t22 - t23 * t6 + t47; t23 * t28 + t24 * t4, -t23 * t27 + t24 * t5, -t52, t23 * t22 + t24 * t6 + t48; t12 - t16 / 0.2e1, t11 + t15 / 0.2e1, t39, pkin(3) * t51 + t39 * t43 + t45; -t23 * t2 + t24 * t20, -t24 * t19 - t23 * t3, t10, -t23 * t1 + t24 * t7 + t47; t24 * t2 + t23 * t20, -t23 * t19 + t24 * t3, -t52, t24 * t1 + t23 * t7 + t48; t9 - t14 / 0.2e1, t8 + t13 / 0.2e1, t39, t39 * t34 + t44 * t37 + t45;];
Tc_stack = t21;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
