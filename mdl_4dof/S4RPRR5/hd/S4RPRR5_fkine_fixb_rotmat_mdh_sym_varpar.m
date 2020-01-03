% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRR5
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RPRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:28
% EndTime: 2019-12-31 16:51:28
% DurationCPUTime: 0.07s
% Computational Cost: add. (46->24), mult. (52->18), div. (0->0), fcn. (86->6), ass. (0->15)
t22 = cos(qJ(3));
t21 = sin(qJ(1));
t20 = sin(qJ(3));
t12 = pkin(4) + 0;
t15 = cos(qJ(1));
t19 = t15 * pkin(1) + t21 * qJ(2) + 0;
t18 = t15 * pkin(2) + t19;
t17 = t21 * pkin(1) - t15 * qJ(2) + 0;
t16 = t21 * pkin(2) + t17;
t14 = cos(qJ(4));
t13 = sin(qJ(4));
t7 = -pkin(5) + t12;
t2 = t15 * t20 - t21 * t22;
t1 = -t15 * t22 - t21 * t20;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t15, -t21, 0, 0; t21, t15, 0, 0; 0, 0, 1, t12; 0, 0, 0, 1; t15, 0, t21, t19; t21, 0, -t15, t17; 0, 1, 0, t12; 0, 0, 0, 1; -t1, -t2, 0, t18; -t2, t1, 0, t16; 0, 0, -1, t7; 0, 0, 0, 1; -t1 * t14, t1 * t13, t2, -t1 * pkin(3) + t2 * pkin(6) + t18; -t2 * t14, t2 * t13, -t1, -t2 * pkin(3) - t1 * pkin(6) + t16; -t13, -t14, 0, t7; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
