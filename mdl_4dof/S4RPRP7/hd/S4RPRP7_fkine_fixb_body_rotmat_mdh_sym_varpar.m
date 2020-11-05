% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:42
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RPRP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:42:50
	% EndTime: 2020-11-04 19:42:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:42:50
	% EndTime: 2020-11-04 19:42:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t23 = cos(qJ(1));
	t22 = sin(qJ(1));
	t1 = [t23, -t22, 0, 0; t22, t23, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:42:50
	% EndTime: 2020-11-04 19:42:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t25 = cos(qJ(1));
	t24 = sin(qJ(1));
	t1 = [0, -t25, t24, t25 * pkin(1) + t24 * qJ(2) + 0; 0, -t24, -t25, t24 * pkin(1) - t25 * qJ(2) + 0; 1, 0, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:42:50
	% EndTime: 2020-11-04 19:42:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t30 = pkin(1) + pkin(5);
	t29 = cos(qJ(1));
	t28 = cos(qJ(3));
	t27 = sin(qJ(1));
	t26 = sin(qJ(3));
	t1 = [t27 * t26, t27 * t28, t29, t27 * qJ(2) + t30 * t29 + 0; -t29 * t26, -t29 * t28, t27, -t29 * qJ(2) + t30 * t27 + 0; t28, -t26, 0, pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:42:50
	% EndTime: 2020-11-04 19:42:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->13), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->7)
	t36 = pkin(1) + pkin(5);
	t35 = cos(qJ(1));
	t34 = cos(qJ(3));
	t33 = sin(qJ(1));
	t32 = sin(qJ(3));
	t31 = -t32 * pkin(3) + t34 * qJ(4) - qJ(2);
	t1 = [t33 * t32, t35, -t33 * t34, -t31 * t33 + t36 * t35 + 0; -t35 * t32, t33, t35 * t34, t31 * t35 + t36 * t33 + 0; t34, 0, t32, t34 * pkin(3) + t32 * qJ(4) + pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end