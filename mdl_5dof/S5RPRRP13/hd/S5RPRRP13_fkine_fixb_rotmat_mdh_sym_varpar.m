% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRRP13_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:25
% EndTime: 2019-12-31 18:58:25
% DurationCPUTime: 0.11s
% Computational Cost: add. (70->38), mult. (98->34), div. (0->0), fcn. (144->6), ass. (0->29)
t19 = sin(qJ(3));
t20 = sin(qJ(1));
t37 = t20 * t19;
t21 = cos(qJ(4));
t36 = t20 * t21;
t22 = cos(qJ(3));
t35 = t20 * t22;
t18 = sin(qJ(4));
t34 = t22 * t18;
t23 = cos(qJ(1));
t33 = t23 * t19;
t32 = t23 * t21;
t9 = t23 * t22;
t17 = pkin(5) + 0;
t31 = t20 * pkin(1) + 0;
t30 = pkin(2) + t17;
t29 = t23 * pkin(1) + t20 * qJ(2) + 0;
t28 = t23 * pkin(6) + t29;
t27 = t22 * pkin(3) + t19 * pkin(7) + t30;
t26 = -t23 * qJ(2) + t31;
t25 = pkin(3) * t37 - pkin(7) * t35 + t28;
t12 = t20 * pkin(6);
t24 = t12 + pkin(7) * t9 + (-pkin(3) * t19 - qJ(2)) * t23 + t31;
t7 = t22 * t21;
t4 = t20 * t18 - t19 * t32;
t3 = t18 * t33 + t36;
t2 = t23 * t18 + t19 * t36;
t1 = t18 * t37 - t32;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t20, 0, 0; t20, t23, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; 0, -t23, t20, t29; 0, -t20, -t23, t26; 1, 0, 0, t17; 0, 0, 0, 1; t37, t35, t23, t28; -t33, -t9, t20, t12 + t26; t22, -t19, 0, t30; 0, 0, 0, 1; t2, -t1, -t35, t25; t4, t3, t9, t24; t7, -t34, t19, t27; 0, 0, 0, 1; t2, -t35, t1, t2 * pkin(4) + t1 * qJ(5) + t25; t4, t9, -t3, t4 * pkin(4) - t3 * qJ(5) + t24; t7, t19, t34, (pkin(4) * t21 + qJ(5) * t18) * t22 + t27; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
