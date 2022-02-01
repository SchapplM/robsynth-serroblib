% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPRPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:45
% EndTime: 2022-01-23 09:24:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (120->58), mult. (110->64), div. (0->0), fcn. (161->10), ass. (0->39)
t24 = sin(qJ(3));
t42 = t24 * pkin(3);
t26 = cos(qJ(3));
t41 = t26 * pkin(3);
t22 = cos(pkin(8));
t40 = pkin(2) * t22 + pkin(1);
t21 = sin(pkin(8));
t39 = t21 * t26;
t25 = sin(qJ(1));
t38 = t25 * t22;
t37 = t25 * t24;
t36 = t25 * t26;
t27 = cos(qJ(1));
t35 = t27 * t22;
t34 = t27 * t24;
t33 = t27 * t26;
t23 = qJ(4) + pkin(6);
t20 = pkin(5) + 0;
t19 = qJ(3) + pkin(9);
t32 = t25 * qJ(2) + 0;
t31 = t21 * pkin(2) + t20;
t30 = t27 * pkin(1) + t32;
t29 = -t27 * qJ(2) + 0;
t18 = -pkin(7) - t23;
t11 = cos(t19);
t2 = pkin(4) * t11 + pkin(2) + t41;
t28 = -t18 * t21 + t2 * t22;
t16 = t25 * pkin(1);
t12 = qJ(5) + t19;
t10 = sin(t19);
t9 = qJ(2) + t42;
t8 = cos(t12);
t7 = sin(t12);
t6 = t27 * t21;
t5 = t25 * t21;
t4 = t21 * pkin(6) + t40;
t3 = pkin(4) * t10 + t42;
t1 = t23 * t21 + t22 * t41 + t40;
t13 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t27, -t25, 0, 0; t25, t27, 0, 0; 0, 0, 1, t20; t35, -t6, t25, t30; t38, -t5, -t27, t16 + t29; t21, t22, 0, t20; t22 * t33 + t37, -t22 * t34 + t36, t6, t4 * t27 + t32; t22 * t36 - t34, -t22 * t37 - t33, t5, t4 * t25 + t29; t39, -t21 * t24, -t22, -t22 * pkin(6) + t31; t25 * t10 + t11 * t35, -t10 * t35 + t25 * t11, t6, t1 * t27 + t9 * t25 + 0; -t27 * t10 + t11 * t38, -t10 * t38 - t27 * t11, t5, t1 * t25 - t9 * t27 + 0; t21 * t11, -t21 * t10, -t22, pkin(3) * t39 - t22 * t23 + t31; t25 * t7 + t8 * t35, t25 * t8 - t7 * t35, t6, t25 * t3 + t27 * t28 + t30; -t27 * t7 + t8 * t38, -t27 * t8 - t7 * t38, t5, t16 + 0 + (-qJ(2) - t3) * t27 + t28 * t25; t21 * t8, -t21 * t7, -t22, t22 * t18 + t21 * t2 + t20;];
Tc_stack = t13;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
