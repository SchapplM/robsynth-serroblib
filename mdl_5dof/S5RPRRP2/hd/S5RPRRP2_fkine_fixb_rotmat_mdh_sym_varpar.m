% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPRRP2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:27:45
% EndTime: 2022-01-23 09:27:45
% DurationCPUTime: 0.09s
% Computational Cost: add. (106->29), mult. (33->18), div. (0->0), fcn. (61->8), ass. (0->24)
t16 = sin(qJ(4));
t14 = qJ(1) + pkin(8);
t11 = qJ(3) + t14;
t5 = sin(t11);
t27 = t5 * t16;
t6 = cos(t11);
t26 = t6 * t16;
t25 = pkin(5) + 0;
t17 = sin(qJ(1));
t24 = t17 * pkin(1) + 0;
t19 = cos(qJ(1));
t23 = t19 * pkin(1) + 0;
t9 = sin(t14);
t22 = pkin(2) * t9 + t24;
t10 = cos(t14);
t21 = pkin(2) * t10 + t23;
t20 = qJ(2) + t25;
t8 = pkin(6) + t20;
t18 = cos(qJ(4));
t15 = -qJ(5) - pkin(7);
t7 = pkin(4) * t18 + pkin(3);
t2 = t6 * t18;
t1 = t5 * t18;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t19, -t17, 0, 0; t17, t19, 0, 0; 0, 0, 1, t25; t10, -t9, 0, t23; t9, t10, 0, t24; 0, 0, 1, t20; t6, -t5, 0, t21; t5, t6, 0, t22; 0, 0, 1, t8; t2, -t26, t5, pkin(3) * t6 + pkin(7) * t5 + t21; t1, -t27, -t6, pkin(3) * t5 - pkin(7) * t6 + t22; t16, t18, 0, t8; t2, -t26, t5, -t15 * t5 + t6 * t7 + t21; t1, -t27, -t6, t15 * t6 + t5 * t7 + t22; t16, t18, 0, pkin(4) * t16 + t8;];
Tc_stack = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
