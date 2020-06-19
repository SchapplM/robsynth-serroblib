% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% Tc_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)
% T_c_stack [(6+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 10:07
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S6RRRRRR10V2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 10:07:09
% EndTime: 2020-06-19 10:07:09
% DurationCPUTime: 0.32s
% Computational Cost: add. (182->47), mult. (202->57), div. (0->0), fcn. (295->12), ass. (0->42)
t28 = qJ(2) + qJ(3);
t24 = sin(t28);
t31 = sin(qJ(4));
t51 = t24 * t31;
t36 = cos(qJ(4));
t50 = t24 * t36;
t33 = sin(qJ(1));
t49 = t33 * t24;
t25 = cos(t28);
t48 = t33 * t25;
t47 = t33 * t31;
t46 = t33 * t36;
t38 = cos(qJ(1));
t45 = t38 * t24;
t44 = t38 * t25;
t43 = t38 * t31;
t42 = t38 * t36;
t27 = pkin(4) + 0;
t37 = cos(qJ(2));
t23 = t37 * pkin(2) + pkin(1);
t41 = t33 * t23 + 0;
t40 = t38 * t23 + 0;
t32 = sin(qJ(2));
t39 = t32 * pkin(2) + t27;
t5 = pkin(3) * t48 + pkin(5) * t49 + t41;
t6 = pkin(3) * t44 + pkin(5) * t45 + t40;
t7 = t24 * pkin(3) - t25 * pkin(5) + t39;
t35 = cos(qJ(5));
t34 = cos(qJ(6));
t30 = sin(qJ(5));
t29 = sin(qJ(6));
t13 = t25 * t42 + t47;
t12 = t25 * t43 - t46;
t11 = t25 * t46 - t43;
t10 = t25 * t47 + t42;
t9 = -t25 * t30 + t35 * t50;
t8 = t25 * t35 + t30 * t50;
t4 = t13 * t35 + t30 * t45;
t3 = t13 * t30 - t35 * t45;
t2 = t11 * t35 + t30 * t49;
t1 = t11 * t30 - t35 * t49;
t14 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t38, -t33, 0, 0; t33, t38, 0, 0; 0, 0, 1, t27; t38 * t37, -t38 * t32, t33, t38 * pkin(1) + 0; t33 * t37, -t33 * t32, -t38, t33 * pkin(1) + 0; t32, t37, 0, t27; t44, -t45, t33, t40; t48, -t49, -t38, t41; t24, t25, 0, t39; t13, -t12, t45, t6; t11, -t10, t49, t5; t50, -t51, -t25, t7; t4, -t3, t12, t6; t2, -t1, t10, t5; t9, -t8, t51, t7; t12 * t29 + t4 * t34, t12 * t34 - t4 * t29, t3, t3 * pkin(6) + t6; t10 * t29 + t2 * t34, t10 * t34 - t2 * t29, t1, t1 * pkin(6) + t5; t29 * t51 + t9 * t34, -t9 * t29 + t34 * t51, t8, t8 * pkin(6) + t7;];
Tc_stack = t14;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,6+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
