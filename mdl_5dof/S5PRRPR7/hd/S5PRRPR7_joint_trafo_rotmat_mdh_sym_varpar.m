% Calculate homogenous joint transformation matrices for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-02-01 03:10
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5PRRPR7_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-02-01 03:10:26
% EndTime: 2022-02-01 03:10:26
% DurationCPUTime: 0.03s
% Computational Cost: add. (10->10), mult. (6->6), div. (0->0), fcn. (30->12), ass. (0->13)
t102 = cos(qJ(2));
t101 = cos(qJ(3));
t100 = cos(qJ(5));
t99 = sin(qJ(2));
t98 = sin(qJ(3));
t97 = sin(qJ(5));
t96 = cos(pkin(5));
t95 = cos(pkin(9));
t94 = cos(pkin(10));
t93 = sin(pkin(5));
t92 = sin(pkin(9));
t91 = sin(pkin(10));
t1 = [t95, -t92, 0, 0; t92, t95, 0, 0; 0, 0, 1, qJ(1); t102, -t99, 0, pkin(1); t96 * t99, t96 * t102, -t93, -t93 * pkin(6); t93 * t99, t93 * t102, t96, t96 * pkin(6); t101, -t98, 0, pkin(2); 0, 0, -1, -pkin(7); t98, t101, 0, 0; t94, -t91, 0, pkin(3); 0, 0, -1, -qJ(4); t91, t94, 0, 0; t100, -t97, 0, pkin(4); 0, 0, -1, -pkin(8); t97, t100, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
