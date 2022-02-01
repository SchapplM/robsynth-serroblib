% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR6
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:41
% EndTime: 2022-01-20 11:16:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (128->42), mult. (83->48), div. (0->0), fcn. (125->10), ass. (0->31)
t17 = cos(pkin(9));
t15 = qJ(1) + qJ(2);
t8 = sin(t15);
t35 = t8 * t17;
t18 = sin(qJ(4));
t34 = t8 * t18;
t10 = cos(t15);
t33 = t10 * t17;
t32 = t17 * t18;
t20 = cos(qJ(4));
t31 = t17 * t20;
t30 = pkin(5) + 0;
t19 = sin(qJ(1));
t29 = t19 * pkin(1) + 0;
t21 = cos(qJ(1));
t28 = t21 * pkin(1) + 0;
t11 = pkin(6) + t30;
t27 = t8 * pkin(2) + t29;
t26 = t10 * pkin(2) + t8 * qJ(3) + t28;
t16 = sin(pkin(9));
t25 = pkin(3) * t17 + pkin(7) * t16;
t22 = -pkin(8) - pkin(7);
t6 = t20 * pkin(4) + pkin(3);
t24 = -t16 * t22 + t17 * t6;
t23 = -t10 * qJ(3) + t27;
t14 = qJ(4) + qJ(5);
t9 = cos(t14);
t7 = sin(t14);
t2 = t10 * t16;
t1 = t8 * t16;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t21, -t19, 0, 0; t19, t21, 0, 0; 0, 0, 1, t30; t10, -t8, 0, t28; t8, t10, 0, t29; 0, 0, 1, t11; t33, -t2, t8, t26; t35, -t1, -t10, t23; t16, t17, 0, t11; t10 * t31 + t34, -t10 * t32 + t8 * t20, t2, t25 * t10 + t26; -t10 * t18 + t8 * t31, -t10 * t20 - t8 * t32, t1, t25 * t8 + t23; t16 * t20, -t16 * t18, -t17, t16 * pkin(3) - t17 * pkin(7) + t11; t9 * t33 + t8 * t7, -t7 * t33 + t8 * t9, t2, pkin(4) * t34 + t24 * t10 + t26; -t10 * t7 + t9 * t35, -t10 * t9 - t7 * t35, t1, t24 * t8 + (-pkin(4) * t18 - qJ(3)) * t10 + t27; t16 * t9, -t16 * t7, -t17, t16 * t6 + t17 * t22 + t11;];
Tc_stack = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
