% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% T_c_mdh [4x4x(2+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   3:  mdh base (link 0) -> mdh frame (3-1), link (3-1)
%   ...
%   2+1:  mdh base (link 0) -> mdh frame (2)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S2RR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:04
% EndTime: 2020-01-03 11:19:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (8->8), mult. (6->6), div. (0->0), fcn. (18->4), ass. (0->5)
t4 = cos(qJ(1));
t3 = cos(qJ(2));
t2 = sin(qJ(1));
t1 = sin(qJ(2));
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t2, t4, 0, 0; 0, 0, 1, 0; t4, -t2, 0, 0; 0, 0, 0, 1; t2 * t3, -t2 * t1, t4, t4 * pkin(1) + 0; -t1, -t3, 0, 0; t4 * t3, -t4 * t1, -t2, -t2 * pkin(1) + 0; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,2+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,2+1]); end % symbolisch
for i = 1:2+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
