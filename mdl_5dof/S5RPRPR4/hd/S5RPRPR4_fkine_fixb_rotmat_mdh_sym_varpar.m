% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR4
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPRPR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:17
% EndTime: 2022-01-23 09:22:18
% DurationCPUTime: 0.10s
% Computational Cost: add. (110->36), mult. (41->30), div. (0->0), fcn. (73->10), ass. (0->23)
t21 = cos(qJ(3));
t4 = t21 * pkin(3) + pkin(2);
t18 = -qJ(4) - pkin(6);
t26 = pkin(5) + 0;
t16 = qJ(3) + pkin(9);
t20 = sin(qJ(1));
t25 = t20 * pkin(1) + 0;
t22 = cos(qJ(1));
t24 = t22 * pkin(1) + 0;
t10 = qJ(2) + t26;
t19 = sin(qJ(3));
t23 = t19 * pkin(3) + t10;
t17 = qJ(1) + pkin(8);
t15 = pkin(7) - t18;
t9 = qJ(5) + t16;
t8 = cos(t17);
t7 = cos(t16);
t6 = sin(t17);
t5 = sin(t16);
t3 = cos(t9);
t2 = sin(t9);
t1 = pkin(4) * t7 + t4;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t22, -t20, 0, 0; t20, t22, 0, 0; 0, 0, 1, t26; t8, -t6, 0, t24; t6, t8, 0, t25; 0, 0, 1, t10; t8 * t21, -t8 * t19, t6, pkin(2) * t8 + pkin(6) * t6 + t24; t6 * t21, -t6 * t19, -t8, pkin(2) * t6 - pkin(6) * t8 + t25; t19, t21, 0, t10; t8 * t7, -t8 * t5, t6, -t18 * t6 + t4 * t8 + t24; t6 * t7, -t6 * t5, -t8, t18 * t8 + t4 * t6 + t25; t5, t7, 0, t23; t8 * t3, -t8 * t2, t6, t1 * t8 + t15 * t6 + t24; t6 * t3, -t6 * t2, -t8, t1 * t6 - t15 * t8 + t25; t2, t3, 0, pkin(4) * t5 + t23;];
Tc_stack = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
