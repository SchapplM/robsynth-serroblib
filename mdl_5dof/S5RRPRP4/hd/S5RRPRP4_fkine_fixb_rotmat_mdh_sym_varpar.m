% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRPRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:31
% EndTime: 2019-12-31 19:52:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (91->29), mult. (42->18), div. (0->0), fcn. (70->6), ass. (0->23)
t15 = cos(qJ(4));
t12 = qJ(1) + qJ(2);
t7 = sin(t12);
t28 = t7 * t15;
t13 = sin(qJ(4));
t8 = cos(t12);
t27 = t8 * t13;
t26 = t8 * t15;
t25 = pkin(5) + 0;
t14 = sin(qJ(1));
t24 = t14 * pkin(1) + 0;
t16 = cos(qJ(1));
t23 = t16 * pkin(1) + 0;
t9 = pkin(6) + t25;
t22 = t7 * pkin(2) + t24;
t21 = pkin(3) + t9;
t20 = t8 * pkin(2) + t7 * qJ(3) + t23;
t19 = t8 * pkin(7) + t20;
t18 = pkin(4) * t13 - qJ(5) * t15;
t17 = -t8 * qJ(3) + t22;
t3 = t7 * pkin(7);
t1 = t7 * t13;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t14, 0, 0; t14, t16, 0, 0; 0, 0, 1, t25; 0, 0, 0, 1; t8, -t7, 0, t23; t7, t8, 0, t24; 0, 0, 1, t9; 0, 0, 0, 1; 0, -t8, t7, t20; 0, -t7, -t8, t17; 1, 0, 0, t9; 0, 0, 0, 1; t1, t28, t8, t19; -t27, -t26, t7, t17 + t3; t15, -t13, 0, t21; 0, 0, 0, 1; t1, t8, -t28, t18 * t7 + t19; -t27, t7, t26, t3 + (-qJ(3) - t18) * t8 + t22; t15, 0, t13, t15 * pkin(4) + t13 * qJ(5) + t21; 0, 0, 0, 1;];
T_ges = t2;
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
