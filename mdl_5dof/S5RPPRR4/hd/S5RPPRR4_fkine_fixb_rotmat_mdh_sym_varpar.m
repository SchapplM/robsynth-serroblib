% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:08
% EndTime: 2022-01-23 09:16:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (120->58), mult. (107->64), div. (0->0), fcn. (158->10), ass. (0->32)
t21 = sin(pkin(9));
t34 = t21 * pkin(3);
t23 = cos(pkin(9));
t10 = t23 * pkin(3) + pkin(2);
t24 = cos(pkin(8));
t26 = sin(qJ(1));
t33 = t26 * t24;
t27 = cos(qJ(1));
t32 = t27 * t24;
t25 = qJ(3) + pkin(6);
t20 = pkin(5) + 0;
t19 = pkin(9) + qJ(4);
t31 = t26 * qJ(2) + 0;
t30 = t27 * pkin(1) + t31;
t29 = -t27 * qJ(2) + 0;
t18 = -pkin(7) - t25;
t12 = cos(t19);
t2 = pkin(4) * t12 + t10;
t22 = sin(pkin(8));
t28 = -t18 * t22 + t2 * t24;
t16 = t26 * pkin(1);
t13 = qJ(5) + t19;
t11 = sin(t19);
t9 = cos(t13);
t8 = sin(t13);
t7 = t27 * t22;
t6 = t26 * t22;
t5 = qJ(2) + t34;
t4 = pkin(2) * t24 + t22 * qJ(3) + pkin(1);
t3 = pkin(4) * t11 + t34;
t1 = t10 * t24 + t25 * t22 + pkin(1);
t14 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t27, -t26, 0, 0; t26, t27, 0, 0; 0, 0, 1, t20; t32, -t7, t26, t30; t33, -t6, -t27, t16 + t29; t22, t24, 0, t20; t26 * t21 + t23 * t32, -t21 * t32 + t26 * t23, t7, t4 * t27 + t31; -t27 * t21 + t23 * t33, -t21 * t33 - t27 * t23, t6, t4 * t26 + t29; t22 * t23, -t22 * t21, -t24, t22 * pkin(2) - t24 * qJ(3) + t20; t26 * t11 + t12 * t32, -t11 * t32 + t26 * t12, t7, t1 * t27 + t5 * t26 + 0; -t27 * t11 + t12 * t33, -t11 * t33 - t27 * t12, t6, t1 * t26 - t5 * t27 + 0; t22 * t12, -t22 * t11, -t24, t22 * t10 - t24 * t25 + t20; t26 * t8 + t9 * t32, t26 * t9 - t8 * t32, t7, t26 * t3 + t27 * t28 + t30; -t27 * t8 + t9 * t33, -t27 * t9 - t8 * t33, t6, t16 + 0 + (-qJ(2) - t3) * t27 + t28 * t26; t22 * t9, -t22 * t8, -t24, t24 * t18 + t22 * t2 + t20;];
Tc_stack = t14;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
