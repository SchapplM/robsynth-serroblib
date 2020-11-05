% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPRP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:31
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:29
	% EndTime: 2020-11-04 19:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:29
	% EndTime: 2020-11-04 19:31:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 0, 1, qJ(1) + 0; 0, -1, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:29
	% EndTime: 2020-11-04 19:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t19 = cos(pkin(5));
	t18 = sin(pkin(5));
	t1 = [t19, -t18, 0, pkin(1) + 0; t18, t19, 0, qJ(1) + 0; 0, 0, 1, qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:29
	% EndTime: 2020-11-04 19:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (11->8), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t22 = pkin(5) + qJ(3);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [t21, -t20, 0, cos(pkin(5)) * pkin(2) + pkin(1) + 0; t20, t21, 0, sin(pkin(5)) * pkin(2) + qJ(1) + 0; 0, 0, 1, pkin(4) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:29
	% EndTime: 2020-11-04 19:31:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t25 = pkin(5) + qJ(3);
	t24 = cos(t25);
	t23 = sin(t25);
	t1 = [t24, 0, t23, t24 * pkin(3) + t23 * qJ(4) + cos(pkin(5)) * pkin(2) + pkin(1) + 0; t23, 0, -t24, t23 * pkin(3) - t24 * qJ(4) + sin(pkin(5)) * pkin(2) + qJ(1) + 0; 0, 1, 0, pkin(4) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end