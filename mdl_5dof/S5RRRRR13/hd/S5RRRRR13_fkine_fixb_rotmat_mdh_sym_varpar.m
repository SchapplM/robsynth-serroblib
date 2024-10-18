% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RRRRR13_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:10
% EndTime: 2024-09-27 17:30:11
% DurationCPUTime: 0.00s
% Computational Cost: add. (165->43), mult. (74->46), div. (0->0), fcn. (106->16), ass. (0->39)
t24 = qJ(1) + qJ(2);
t20 = qJ(3) + t24;
t11 = sin(t20);
t25 = sin(pkin(5));
t4 = t11 * t25;
t12 = cos(t20);
t42 = t12 * t25;
t27 = sin(qJ(4));
t41 = t25 * t27;
t26 = cos(pkin(5));
t40 = t26 * t27;
t29 = cos(qJ(4));
t39 = t26 * t29;
t23 = qJ(4) + qJ(5);
t38 = pkin(6) + 0;
t28 = sin(qJ(1));
t37 = t28 * pkin(1) + 0;
t30 = cos(qJ(1));
t36 = t30 * pkin(1) + 0;
t35 = pkin(7) + t38;
t17 = sin(t24);
t34 = pkin(2) * t17 + t37;
t19 = cos(t24);
t33 = pkin(2) * t19 + t36;
t32 = pkin(8) + t35;
t31 = pkin(9) + pkin(10);
t18 = cos(t23);
t16 = sin(t23);
t15 = pkin(5) - t23;
t14 = pkin(5) + t23;
t13 = t29 * pkin(4) + pkin(3);
t8 = cos(t14);
t7 = sin(t15);
t6 = cos(t15) / 0.2e1;
t5 = sin(t14) / 0.2e1;
t3 = pkin(4) * t40 - t25 * t31;
t2 = t6 + t8 / 0.2e1;
t1 = t5 - t7 / 0.2e1;
t9 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t30, -t28, 0, 0; t28, t30, 0, 0; 0, 0, 1, t38; t19, -t17, 0, t36; t17, t19, 0, t37; 0, 0, 1, t35; t12, -t11, 0, t33; t11, t12, 0, t34; 0, 0, 1, t32; -t11 * t40 + t12 * t29, -t11 * t39 - t12 * t27, t4, t12 * pkin(3) + pkin(9) * t4 + t33; t11 * t29 + t12 * t40, -t11 * t27 + t12 * t39, -t42, t11 * pkin(3) - pkin(9) * t42 + t34; t41, t25 * t29, t26, t26 * pkin(9) + t32; -t11 * t1 + t12 * t18, -t11 * t2 - t12 * t16, t4, -t11 * t3 + t12 * t13 + t33; t12 * t1 + t11 * t18, -t11 * t16 + t12 * t2, -t42, t11 * t13 + t12 * t3 + t34; t6 - t8 / 0.2e1, t5 + t7 / 0.2e1, t26, pkin(4) * t41 + t26 * t31 + t32;];
Tc_stack = t9;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
