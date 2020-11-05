% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:41
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:41:13
	% EndTime: 2020-11-04 19:41:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:41:13
	% EndTime: 2020-11-04 19:41:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t23 = cos(qJ(1));
	t22 = sin(qJ(1));
	t1 = [t23, -t22, 0, 0; t22, t23, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:41:13
	% EndTime: 2020-11-04 19:41:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t25 = cos(qJ(1));
	t24 = sin(qJ(1));
	t1 = [0, -t25, t24, t25 * pkin(1) + t24 * qJ(2) + 0; 0, -t24, -t25, t24 * pkin(1) - t25 * qJ(2) + 0; 1, 0, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:41:13
	% EndTime: 2020-11-04 19:41:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t30 = cos(qJ(1));
	t29 = sin(qJ(1));
	t28 = pkin(1) + qJ(3);
	t27 = cos(pkin(6));
	t26 = sin(pkin(6));
	t1 = [t29 * t26, t29 * t27, t30, t29 * qJ(2) + t28 * t30 + 0; -t30 * t26, -t30 * t27, t29, -t30 * qJ(2) + t28 * t29 + 0; t27, -t26, 0, pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:41:13
	% EndTime: 2020-11-04 19:41:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t35 = pkin(6) + qJ(4);
	t34 = pkin(1) + pkin(5) + qJ(3);
	t33 = cos(t35);
	t32 = sin(t35);
	t31 = sin(pkin(6)) * pkin(3) + qJ(2);
	t1 = [t36 * t32, t36 * t33, t37, t31 * t36 + t34 * t37 + 0; -t37 * t32, -t37 * t33, t36, -t31 * t37 + t34 * t36 + 0; t33, -t32, 0, cos(pkin(6)) * pkin(3) + pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end