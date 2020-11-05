% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:37
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:58
	% EndTime: 2020-11-04 19:37:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:58
	% EndTime: 2020-11-04 19:37:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, -1, 0, 0; 0, 0, -1, -qJ(1) + 0; 1, 0, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:58
	% EndTime: 2020-11-04 19:37:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t20 = cos(qJ(2));
	t19 = sin(qJ(2));
	t1 = [-t19, -t20, 0, 0; 0, 0, -1, -qJ(1) + 0; t20, -t19, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:58
	% EndTime: 2020-11-04 19:37:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->7), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t23 = qJ(2) + qJ(3);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [-t21, -t22, 0, -sin(qJ(2)) * pkin(1) + 0; 0, 0, -1, -qJ(1) + 0; t22, -t21, 0, cos(qJ(2)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:58
	% EndTime: 2020-11-04 19:37:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t27 = qJ(2) + qJ(3);
	t26 = qJ(4) + t27;
	t25 = cos(t26);
	t24 = sin(t26);
	t1 = [-t24, -t25, 0, -pkin(2) * sin(t27) - sin(qJ(2)) * pkin(1) + 0; 0, 0, -1, -qJ(1) + 0; t25, -t24, 0, pkin(2) * cos(t27) + cos(qJ(2)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end