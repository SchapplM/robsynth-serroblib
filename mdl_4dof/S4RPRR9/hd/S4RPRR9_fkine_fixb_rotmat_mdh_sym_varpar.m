% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-29 13:21
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RPRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 13:19:40
% EndTime: 2019-12-29 13:19:40
% DurationCPUTime: 0.15s
% Computational Cost: add. (41->28), mult. (46->26), div. (0->0), fcn. (75->6), ass. (0->21)
t8 = sin(qJ(3));
t9 = sin(qJ(1));
t24 = t9 * t8;
t12 = cos(qJ(1));
t23 = t12 * t8;
t10 = cos(qJ(4));
t22 = t9 * t10;
t11 = cos(qJ(3));
t21 = t9 * t11;
t20 = t12 * t10;
t19 = t12 * t11;
t6 = pkin(4) + 0;
t18 = t9 * pkin(1) + 0;
t17 = pkin(2) + t6;
t16 = t12 * pkin(1) + t9 * qJ(2) + 0;
t15 = t12 * pkin(5) + t16;
t14 = pkin(3) * t8 - pkin(6) * t11;
t13 = -t12 * qJ(2) + t18;
t7 = sin(qJ(4));
t2 = t9 * pkin(5);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t12, -t9, 0, 0; t9, t12, 0, 0; 0, 0, 1, t6; 0, 0, 0, 1; 0, -t12, t9, t16; 0, -t9, -t12, t13; 1, 0, 0, t6; 0, 0, 0, 1; t24, t21, t12, t15; -t23, -t19, t9, t13 + t2; t11, -t8, 0, t17; 0, 0, 0, 1; t12 * t7 + t8 * t22, -t7 * t24 + t20, -t21, t14 * t9 + t15; -t8 * t20 + t9 * t7, t7 * t23 + t22, t19, t2 + (-qJ(2) - t14) * t12 + t18; t11 * t10, -t11 * t7, t8, t11 * pkin(3) + t8 * pkin(6) + t17; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
