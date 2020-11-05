% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:30
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:30:41
	% EndTime: 2020-11-04 19:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:30:41
	% EndTime: 2020-11-04 19:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 0, 1, qJ(1) + 0; 0, -1, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:30:41
	% EndTime: 2020-11-04 19:30:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t19 = cos(pkin(5));
	t18 = sin(pkin(5));
	t1 = [t19, -t18, 0, pkin(1) + 0; t18, t19, 0, qJ(1) + 0; 0, 0, 1, qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:30:41
	% EndTime: 2020-11-04 19:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t21 = cos(pkin(5));
	t20 = sin(pkin(5));
	t1 = [t21, 0, t20, t21 * pkin(2) + t20 * qJ(3) + pkin(1) + 0; t20, 0, -t21, t20 * pkin(2) - t21 * qJ(3) + qJ(1) + 0; 0, 1, 0, qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:30:41
	% EndTime: 2020-11-04 19:30:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->12), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t28 = pkin(2) + pkin(3);
	t27 = cos(qJ(4));
	t26 = sin(qJ(4));
	t25 = cos(pkin(5));
	t24 = sin(pkin(5));
	t23 = t24 * t27 - t25 * t26;
	t22 = -t24 * t26 - t25 * t27;
	t1 = [-t22, t23, 0, t24 * qJ(3) + t28 * t25 + pkin(1) + 0; t23, t22, 0, -t25 * qJ(3) + t28 * t24 + qJ(1) + 0; 0, 0, -1, -pkin(4) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end