% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4RPRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:41
% EndTime: 2019-12-31 16:44:41
% DurationCPUTime: 0.07s
% Computational Cost: add. (59->27), mult. (40->24), div. (0->0), fcn. (68->6), ass. (0->19)
t15 = sin(qJ(1));
t10 = pkin(6) + qJ(3);
t7 = sin(t10);
t22 = t15 * t7;
t16 = cos(qJ(1));
t21 = t16 * t7;
t11 = pkin(4) + 0;
t14 = -pkin(5) - qJ(2);
t13 = cos(pkin(6));
t5 = t13 * pkin(2) + pkin(1);
t20 = t16 * t14 + t15 * t5 + 0;
t12 = sin(pkin(6));
t19 = t12 * pkin(2) + t11;
t8 = cos(t10);
t18 = pkin(3) * t8 + qJ(4) * t7;
t17 = -t15 * t14 + t16 * t5 + 0;
t4 = t16 * t8;
t3 = t15 * t8;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t15, 0, 0; t15, t16, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t16 * t13, -t16 * t12, t15, t16 * pkin(1) + t15 * qJ(2) + 0; t15 * t13, -t15 * t12, -t16, t15 * pkin(1) - t16 * qJ(2) + 0; t12, t13, 0, t11; 0, 0, 0, 1; t4, -t21, t15, t17; t3, -t22, -t16, t20; t7, t8, 0, t19; 0, 0, 0, 1; t4, t15, t21, t18 * t16 + t17; t3, -t16, t22, t18 * t15 + t20; t7, 0, -t8, t7 * pkin(3) - t8 * qJ(4) + t19; 0, 0, 0, 1;];
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
