% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S3RPP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:29
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S3RPP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),uint8(0),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S3RPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [3x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:29:30
	% EndTime: 2020-11-04 19:29:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:29:30
	% EndTime: 2020-11-04 19:29:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t11 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [t11, -t10, 0, 0; t10, t11, 0, 0; 0, 0, 1, pkin(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:29:30
	% EndTime: 2020-11-04 19:29:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t13 = cos(qJ(1));
	t12 = sin(qJ(1));
	t1 = [0, -t13, t12, t13 * pkin(1) + t12 * qJ(2) + 0; 0, -t12, -t13, t12 * pkin(1) - t13 * qJ(2) + 0; 1, 0, 0, pkin(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:29:30
	% EndTime: 2020-11-04 19:29:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = cos(qJ(1));
	t15 = sin(qJ(1));
	t14 = pkin(1) + qJ(3);
	t1 = [0, t15, t16, qJ(2) * t15 + t14 * t16 + 0; 0, -t16, t15, -qJ(2) * t16 + t14 * t15 + 0; 1, 0, 0, pkin(2) + pkin(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end