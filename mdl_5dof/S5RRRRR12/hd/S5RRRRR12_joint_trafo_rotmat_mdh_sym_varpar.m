% Calculate homogenous joint transformation matrices for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-02-04 03:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5RRRRR12_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-02-04 03:34:44
% EndTime: 2022-02-04 03:34:44
% DurationCPUTime: 0.03s
% Computational Cost: add. (11->11), mult. (12->12), div. (0->0), fcn. (40->14), ass. (0->15)
t135 = cos(qJ(1));
t134 = cos(qJ(2));
t133 = cos(qJ(3));
t132 = cos(qJ(4));
t131 = cos(qJ(5));
t130 = sin(qJ(1));
t129 = sin(qJ(2));
t128 = sin(qJ(3));
t127 = sin(qJ(4));
t126 = sin(qJ(5));
t125 = cos(pkin(5));
t124 = cos(pkin(6));
t123 = sin(pkin(5));
t122 = sin(pkin(6));
t1 = [t135, -t130, 0, 0; t130, t135, 0, 0; 0, 0, 1, pkin(7); t134, -t129, 0, pkin(1); t125 * t129, t125 * t134, -t123, -t123 * pkin(8); t123 * t129, t123 * t134, t125, t125 * pkin(8); t133, -t128, 0, pkin(2); t124 * t128, t124 * t133, -t122, -t122 * pkin(9); t122 * t128, t122 * t133, t124, t124 * pkin(9); t132, -t127, 0, pkin(3); 0, 0, -1, -pkin(10); t127, t132, 0, 0; t131, -t126, 0, pkin(4); 0, 0, -1, -pkin(11); t126, t131, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
