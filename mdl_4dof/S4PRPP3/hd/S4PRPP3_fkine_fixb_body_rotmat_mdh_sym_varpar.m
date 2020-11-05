% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRPP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:56
	% EndTime: 2020-11-04 19:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:56
	% EndTime: 2020-11-04 19:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 0, 1, qJ(1) + 0; 0, -1, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:56
	% EndTime: 2020-11-04 19:33:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = cos(qJ(2));
	t12 = sin(qJ(2));
	t1 = [t13, -t12, 0, pkin(1) + 0; t12, t13, 0, qJ(1) + 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:57
	% EndTime: 2020-11-04 19:33:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t15 = cos(qJ(2));
	t14 = sin(qJ(2));
	t1 = [t15, 0, t14, t15 * pkin(2) + t14 * qJ(3) + pkin(1) + 0; t14, 0, -t15, t14 * pkin(2) - t15 * qJ(3) + qJ(1) + 0; 0, 1, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:57
	% EndTime: 2020-11-04 19:33:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (11->10), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t18 = pkin(2) + pkin(3);
	t17 = cos(qJ(2));
	t16 = sin(qJ(2));
	t1 = [t17, t16, 0, t16 * qJ(3) + t18 * t17 + pkin(1) + 0; t16, -t17, 0, -t17 * qJ(3) + t18 * t16 + qJ(1) + 0; 0, 0, -1, -qJ(4) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end