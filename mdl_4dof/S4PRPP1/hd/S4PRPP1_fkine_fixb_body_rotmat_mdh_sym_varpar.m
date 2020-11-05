% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRPP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRPP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:24
	% EndTime: 2020-11-04 19:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:24
	% EndTime: 2020-11-04 19:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t21 = cos(pkin(5));
	t20 = sin(pkin(5));
	t1 = [t21, -t20, 0, 0; t20, t21, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:24
	% EndTime: 2020-11-04 19:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t24 = pkin(5) + qJ(2);
	t23 = cos(t24);
	t22 = sin(t24);
	t1 = [t23, -t22, 0, cos(pkin(5)) * pkin(1) + 0; t22, t23, 0, sin(pkin(5)) * pkin(1) + 0; 0, 0, 1, pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:24
	% EndTime: 2020-11-04 19:33:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t27 = pkin(5) + qJ(2);
	t26 = cos(t27);
	t25 = sin(t27);
	t1 = [0, -t26, t25, t26 * pkin(2) + t25 * qJ(3) + cos(pkin(5)) * pkin(1) + 0; 0, -t25, -t26, t25 * pkin(2) - t26 * qJ(3) + sin(pkin(5)) * pkin(1) + 0; 1, 0, 0, pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:24
	% EndTime: 2020-11-04 19:33:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (20->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t31 = pkin(2) + qJ(4);
	t30 = pkin(5) + qJ(2);
	t29 = cos(t30);
	t28 = sin(t30);
	t1 = [0, t28, t29, t31 * t29 + t28 * qJ(3) + cos(pkin(5)) * pkin(1) + 0; 0, -t29, t28, t31 * t28 - t29 * qJ(3) + sin(pkin(5)) * pkin(1) + 0; 1, 0, 0, pkin(3) + pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end