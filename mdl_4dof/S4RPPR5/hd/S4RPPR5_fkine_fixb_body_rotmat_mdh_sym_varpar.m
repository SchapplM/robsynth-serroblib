% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPPR5 (for one body)
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
% Datum: 2020-11-04 19:40
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:40:45
	% EndTime: 2020-11-04 19:40:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:40:45
	% EndTime: 2020-11-04 19:40:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [t27, -t26, 0, 0; t26, t27, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:40:45
	% EndTime: 2020-11-04 19:40:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t29 = cos(qJ(1));
	t28 = sin(qJ(1));
	t1 = [t29, 0, t28, t29 * pkin(1) + t28 * qJ(2) + 0; t28, 0, -t29, t28 * pkin(1) - t29 * qJ(2) + 0; 0, 1, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:40:45
	% EndTime: 2020-11-04 19:40:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t36 = pkin(1) + pkin(2);
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t33 = cos(pkin(6));
	t32 = sin(pkin(6));
	t31 = -t35 * t32 + t34 * t33;
	t30 = -t34 * t32 - t35 * t33;
	t1 = [-t30, t31, 0, t34 * qJ(2) + t36 * t35 + 0; t31, t30, 0, -t35 * qJ(2) + t36 * t34 + 0; 0, 0, -1, -qJ(3) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:40:45
	% EndTime: 2020-11-04 19:40:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (26->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t47 = cos(qJ(1));
	t46 = cos(qJ(4));
	t45 = sin(qJ(1));
	t44 = sin(qJ(4));
	t43 = cos(pkin(6));
	t42 = sin(pkin(6));
	t40 = -pkin(3) * t42 + pkin(5) * t43 - qJ(2);
	t39 = pkin(3) * t43 + pkin(5) * t42 + pkin(1) + pkin(2);
	t38 = t42 * t45 + t43 * t47;
	t37 = t42 * t47 - t43 * t45;
	t1 = [t38 * t46, -t38 * t44, t37, t39 * t47 - t40 * t45 + 0; -t37 * t46, t37 * t44, t38, t39 * t45 + t40 * t47 + 0; -t44, -t46, 0, -qJ(3) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end