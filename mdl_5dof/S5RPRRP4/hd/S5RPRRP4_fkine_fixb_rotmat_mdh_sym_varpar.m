% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP4
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:39
% EndTime: 2022-01-23 09:31:40
% DurationCPUTime: 0.11s
% Computational Cost: add. (110->53), mult. (110->54), div. (0->0), fcn. (161->8), ass. (0->42)
t29 = pkin(7) + pkin(6);
t25 = sin(qJ(3));
t45 = t25 * pkin(3);
t27 = cos(qJ(3));
t44 = t27 * pkin(3);
t24 = cos(pkin(8));
t43 = pkin(2) * t24 + pkin(1);
t22 = qJ(3) + qJ(4);
t13 = sin(t22);
t23 = sin(pkin(8));
t42 = t23 * t13;
t41 = t23 * t27;
t26 = sin(qJ(1));
t40 = t26 * t24;
t39 = t26 * t25;
t38 = t26 * t27;
t28 = cos(qJ(1));
t37 = t28 * t24;
t36 = t28 * t25;
t35 = t28 * t27;
t21 = pkin(5) + 0;
t34 = t26 * qJ(2) + 0;
t33 = t23 * pkin(2) + t21;
t32 = t28 * pkin(1) + t34;
t31 = -t28 * qJ(2) + 0;
t20 = -qJ(5) - t29;
t14 = cos(t22);
t6 = pkin(4) * t14 + pkin(2) + t44;
t30 = -t20 * t23 + t24 * t6;
t18 = t26 * pkin(1);
t12 = qJ(2) + t45;
t11 = t28 * t23;
t10 = t26 * t23;
t9 = t23 * t14;
t8 = t23 * pkin(6) + t43;
t7 = pkin(4) * t13 + t45;
t5 = t26 * t13 + t14 * t37;
t4 = -t13 * t37 + t26 * t14;
t3 = -t28 * t13 + t14 * t40;
t2 = -t13 * t40 - t28 * t14;
t1 = t29 * t23 + t24 * t44 + t43;
t15 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t28, -t26, 0, 0; t26, t28, 0, 0; 0, 0, 1, t21; t37, -t11, t26, t32; t40, -t10, -t28, t18 + t31; t23, t24, 0, t21; t24 * t35 + t39, -t24 * t36 + t38, t11, t8 * t28 + t34; t24 * t38 - t36, -t24 * t39 - t35, t10, t8 * t26 + t31; t41, -t23 * t25, -t24, -t24 * pkin(6) + t33; t5, t4, t11, t1 * t28 + t12 * t26 + 0; t3, t2, t10, t1 * t26 - t12 * t28 + 0; t9, -t42, -t24, pkin(3) * t41 - t24 * t29 + t33; t5, t4, t11, t26 * t7 + t30 * t28 + t32; t3, t2, t10, t18 + 0 + (-qJ(2) - t7) * t28 + t30 * t26; t9, -t42, -t24, t24 * t20 + t23 * t6 + t21;];
Tc_stack = t15;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
