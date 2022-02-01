% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPRRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:22
% EndTime: 2022-01-23 09:34:22
% DurationCPUTime: 0.09s
% Computational Cost: add. (111->24), mult. (26->14), div. (0->0), fcn. (50->10), ass. (0->24)
t29 = pkin(5) + 0;
t16 = qJ(1) + pkin(9);
t18 = sin(qJ(1));
t28 = t18 * pkin(1) + 0;
t20 = cos(qJ(1));
t27 = t20 * pkin(1) + 0;
t11 = sin(t16);
t26 = pkin(2) * t11 + t28;
t12 = cos(t16);
t25 = pkin(2) * t12 + t27;
t13 = qJ(3) + t16;
t24 = qJ(2) + t29;
t8 = sin(t13);
t23 = pkin(3) * t8 + t26;
t9 = cos(t13);
t22 = pkin(3) * t9 + t25;
t21 = pkin(6) + t24;
t19 = cos(qJ(5));
t17 = sin(qJ(5));
t10 = qJ(4) + t13;
t5 = pkin(7) + t21;
t4 = cos(t10);
t3 = sin(t10);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t20, -t18, 0, 0; t18, t20, 0, 0; 0, 0, 1, t29; t12, -t11, 0, t27; t11, t12, 0, t28; 0, 0, 1, t24; t9, -t8, 0, t25; t8, t9, 0, t26; 0, 0, 1, t21; t4, -t3, 0, t22; t3, t4, 0, t23; 0, 0, 1, t5; t4 * t19, -t4 * t17, t3, pkin(4) * t4 + pkin(8) * t3 + t22; t3 * t19, -t3 * t17, -t4, pkin(4) * t3 - pkin(8) * t4 + t23; t17, t19, 0, t5;];
Tc_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
