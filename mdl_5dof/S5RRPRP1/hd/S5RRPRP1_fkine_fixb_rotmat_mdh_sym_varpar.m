% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RRPRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:26
% EndTime: 2022-01-20 10:19:26
% DurationCPUTime: 0.10s
% Computational Cost: add. (106->29), mult. (33->18), div. (0->0), fcn. (61->8), ass. (0->24)
t16 = sin(qJ(4));
t14 = qJ(1) + qJ(2);
t9 = pkin(8) + t14;
t3 = sin(t9);
t27 = t3 * t16;
t4 = cos(t9);
t26 = t4 * t16;
t25 = pkin(5) + 0;
t17 = sin(qJ(1));
t24 = t17 * pkin(1) + 0;
t19 = cos(qJ(1));
t23 = t19 * pkin(1) + 0;
t22 = pkin(6) + t25;
t10 = sin(t14);
t21 = pkin(2) * t10 + t24;
t11 = cos(t14);
t20 = pkin(2) * t11 + t23;
t8 = qJ(3) + t22;
t18 = cos(qJ(4));
t15 = -qJ(5) - pkin(7);
t7 = pkin(4) * t18 + pkin(3);
t2 = t4 * t18;
t1 = t3 * t18;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t19, -t17, 0, 0; t17, t19, 0, 0; 0, 0, 1, t25; t11, -t10, 0, t23; t10, t11, 0, t24; 0, 0, 1, t22; t4, -t3, 0, t20; t3, t4, 0, t21; 0, 0, 1, t8; t2, -t26, t3, pkin(3) * t4 + pkin(7) * t3 + t20; t1, -t27, -t4, pkin(3) * t3 - pkin(7) * t4 + t21; t16, t18, 0, t8; t2, -t26, t3, -t15 * t3 + t4 * t7 + t20; t1, -t27, -t4, t15 * t4 + t3 * t7 + t21; t16, t18, 0, pkin(4) * t16 + t8;];
Tc_stack = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
