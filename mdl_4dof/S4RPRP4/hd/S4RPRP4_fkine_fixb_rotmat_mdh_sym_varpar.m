% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRP4
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
% Datum: 2019-12-29 12:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RPRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 12:51:48
% EndTime: 2019-12-29 12:51:48
% DurationCPUTime: 0.14s
% Computational Cost: add. (62->21), mult. (32->16), div. (0->0), fcn. (56->6), ass. (0->19)
t12 = sin(qJ(3));
t11 = qJ(1) + pkin(6);
t6 = sin(t11);
t23 = t6 * t12;
t7 = cos(t11);
t22 = t7 * t12;
t21 = pkin(4) + 0;
t13 = sin(qJ(1));
t20 = t13 * pkin(1) + 0;
t15 = cos(qJ(1));
t19 = t15 * pkin(1) + 0;
t8 = qJ(2) + t21;
t18 = t7 * pkin(2) + t6 * pkin(5) + t19;
t14 = cos(qJ(3));
t17 = pkin(3) * t14 + qJ(4) * t12;
t16 = t6 * pkin(2) - t7 * pkin(5) + t20;
t2 = t7 * t14;
t1 = t6 * t14;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t15, -t13, 0, 0; t13, t15, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t7, -t6, 0, t19; t6, t7, 0, t20; 0, 0, 1, t8; 0, 0, 0, 1; t2, -t22, t6, t18; t1, -t23, -t7, t16; t12, t14, 0, t8; 0, 0, 0, 1; t2, t6, t22, t17 * t7 + t18; t1, -t7, t23, t17 * t6 + t16; t12, 0, -t14, t12 * pkin(3) - t14 * qJ(4) + t8; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
