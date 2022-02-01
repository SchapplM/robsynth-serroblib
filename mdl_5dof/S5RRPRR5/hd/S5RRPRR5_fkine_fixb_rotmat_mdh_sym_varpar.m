% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RRPRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:04
% EndTime: 2022-01-20 11:02:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (110->36), mult. (41->30), div. (0->0), fcn. (73->10), ass. (0->23)
t19 = cos(pkin(9));
t4 = t19 * pkin(3) + pkin(2);
t20 = -pkin(7) - qJ(3);
t26 = pkin(5) + 0;
t16 = pkin(9) + qJ(4);
t21 = sin(qJ(1));
t25 = t21 * pkin(1) + 0;
t22 = cos(qJ(1));
t24 = t22 * pkin(1) + 0;
t12 = pkin(6) + t26;
t18 = sin(pkin(9));
t23 = t18 * pkin(3) + t12;
t17 = qJ(1) + qJ(2);
t15 = pkin(8) - t20;
t9 = cos(t17);
t8 = sin(t17);
t7 = qJ(5) + t16;
t6 = cos(t16);
t5 = sin(t16);
t3 = cos(t7);
t2 = sin(t7);
t1 = pkin(4) * t6 + t4;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t22, -t21, 0, 0; t21, t22, 0, 0; 0, 0, 1, t26; t9, -t8, 0, t24; t8, t9, 0, t25; 0, 0, 1, t12; t9 * t19, -t9 * t18, t8, t9 * pkin(2) + t8 * qJ(3) + t24; t8 * t19, -t8 * t18, -t9, t8 * pkin(2) - t9 * qJ(3) + t25; t18, t19, 0, t12; t9 * t6, -t9 * t5, t8, -t8 * t20 + t9 * t4 + t24; t8 * t6, -t8 * t5, -t9, t9 * t20 + t8 * t4 + t25; t5, t6, 0, t23; t9 * t3, -t9 * t2, t8, t9 * t1 + t15 * t8 + t24; t8 * t3, -t8 * t2, -t9, t8 * t1 - t9 * t15 + t25; t2, t3, 0, pkin(4) * t5 + t23;];
Tc_stack = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
