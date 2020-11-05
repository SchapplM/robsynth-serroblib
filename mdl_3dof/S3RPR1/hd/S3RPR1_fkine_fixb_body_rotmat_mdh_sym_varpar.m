% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S3RPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:29
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S3RPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S3RPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:29:42
	% EndTime: 2020-11-04 19:29:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:29:42
	% EndTime: 2020-11-04 19:29:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t17 = cos(qJ(1));
	t16 = sin(qJ(1));
	t1 = [t17, -t16, 0, 0; t16, t17, 0, 0; 0, 0, 1, pkin(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:29:42
	% EndTime: 2020-11-04 19:29:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t19 = cos(qJ(1));
	t18 = sin(qJ(1));
	t1 = [t19, 0, t18, pkin(1) * t19 + qJ(2) * t18 + 0; t18, 0, -t19, pkin(1) * t18 - qJ(2) * t19 + 0; 0, 1, 0, pkin(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:29:42
	% EndTime: 2020-11-04 19:29:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t26 = pkin(1) + pkin(2);
	t25 = cos(qJ(1));
	t24 = cos(qJ(3));
	t23 = sin(qJ(1));
	t22 = sin(qJ(3));
	t21 = -t25 * t22 + t23 * t24;
	t20 = -t23 * t22 - t25 * t24;
	t1 = [-t20, t21, 0, t23 * qJ(2) + t26 * t25 + 0; t21, t20, 0, -t25 * qJ(2) + t26 * t23 + 0; 0, 0, -1, -pkin(4) + pkin(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end