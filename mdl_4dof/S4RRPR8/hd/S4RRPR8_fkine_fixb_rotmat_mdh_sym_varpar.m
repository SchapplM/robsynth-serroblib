% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4RRPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:42
% EndTime: 2019-12-31 17:07:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (46->28), mult. (67->28), div. (0->0), fcn. (101->6), ass. (0->21)
t14 = sin(qJ(1));
t13 = sin(qJ(2));
t25 = qJ(3) * t13;
t16 = cos(qJ(2));
t5 = t14 * t16;
t28 = pkin(2) * t5 + t14 * t25;
t27 = t14 * t13;
t17 = cos(qJ(1));
t26 = t17 * t13;
t6 = t17 * t16;
t11 = pkin(4) + 0;
t24 = t14 * pkin(1) + 0;
t23 = t17 * pkin(1) + t14 * pkin(5) + 0;
t22 = -t17 * pkin(5) + t24;
t21 = pkin(2) * t6 + t17 * t25 + t23;
t12 = sin(qJ(4));
t15 = cos(qJ(4));
t20 = -t16 * t12 + t13 * t15;
t19 = t13 * t12 + t16 * t15;
t18 = t13 * pkin(2) - t16 * qJ(3) + t11;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t14, 0, 0; t14, t17, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t6, -t26, t14, t23; t5, -t27, -t17, t22; t13, t16, 0, t11; 0, 0, 0, 1; t6, t14, t26, t21; t5, -t17, t27, t22 + t28; t13, 0, -t16, t18; 0, 0, 0, 1; t19 * t17, t20 * t17, -t14, pkin(3) * t6 - t14 * pkin(6) + t21; t19 * t14, t20 * t14, t17, pkin(3) * t5 + (-pkin(5) + pkin(6)) * t17 + t24 + t28; t20, -t19, 0, t13 * pkin(3) + t18; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
