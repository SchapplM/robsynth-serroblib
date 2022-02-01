% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RRPPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:12
% EndTime: 2022-01-20 09:51:12
% DurationCPUTime: 0.11s
% Computational Cost: add. (112->30), mult. (33->22), div. (0->0), fcn. (61->10), ass. (0->23)
t15 = qJ(1) + qJ(2);
t26 = pkin(5) + 0;
t19 = sin(qJ(1));
t25 = t19 * pkin(1) + 0;
t20 = cos(qJ(1));
t24 = t20 * pkin(1) + 0;
t23 = pkin(6) + t26;
t10 = sin(t15);
t22 = pkin(2) * t10 + t25;
t11 = cos(t15);
t21 = pkin(2) * t11 + t24;
t6 = qJ(3) + t23;
t18 = -pkin(7) - qJ(4);
t17 = cos(pkin(9));
t16 = sin(pkin(9));
t14 = pkin(9) + qJ(5);
t9 = pkin(8) + t15;
t8 = cos(t14);
t7 = sin(t14);
t3 = pkin(4) * t17 + pkin(3);
t2 = cos(t9);
t1 = sin(t9);
t4 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t20, -t19, 0, 0; t19, t20, 0, 0; 0, 0, 1, t26; t11, -t10, 0, t24; t10, t11, 0, t25; 0, 0, 1, t23; t2, -t1, 0, t21; t1, t2, 0, t22; 0, 0, 1, t6; t2 * t17, -t2 * t16, t1, pkin(3) * t2 + qJ(4) * t1 + t21; t1 * t17, -t1 * t16, -t2, pkin(3) * t1 - qJ(4) * t2 + t22; t16, t17, 0, t6; t2 * t8, -t2 * t7, t1, -t1 * t18 + t2 * t3 + t21; t1 * t8, -t1 * t7, -t2, t1 * t3 + t18 * t2 + t22; t7, t8, 0, pkin(4) * t16 + t6;];
Tc_stack = t4;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
