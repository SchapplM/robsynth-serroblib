% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:18
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:59
	% EndTime: 2020-11-04 20:17:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:59
	% EndTime: 2020-11-04 20:17:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t29 = cos(qJ(1));
	t28 = sin(qJ(1));
	t1 = [t29, -t28, 0, 0; t28, t29, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:59
	% EndTime: 2020-11-04 20:17:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t31 = cos(qJ(1));
	t30 = sin(qJ(1));
	t1 = [0, -t31, t30, t31 * pkin(1) + t30 * qJ(2) + 0; 0, -t30, -t31, t30 * pkin(1) - t31 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:59
	% EndTime: 2020-11-04 20:17:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t36 = pkin(1) + pkin(6);
	t35 = cos(qJ(1));
	t34 = cos(qJ(3));
	t33 = sin(qJ(1));
	t32 = sin(qJ(3));
	t1 = [t33 * t32, t33 * t34, t35, t33 * qJ(2) + t36 * t35 + 0; -t35 * t32, -t35 * t34, t33, -t35 * qJ(2) + t36 * t33 + 0; t34, -t32, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:59
	% EndTime: 2020-11-04 20:17:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->13), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->7)
	t42 = pkin(1) + pkin(6);
	t41 = cos(qJ(1));
	t40 = cos(qJ(3));
	t39 = sin(qJ(1));
	t38 = sin(qJ(3));
	t37 = -t38 * pkin(3) + t40 * qJ(4) - qJ(2);
	t1 = [t39 * t38, t41, -t39 * t40, -t37 * t39 + t42 * t41 + 0; -t41 * t38, t39, t41 * t40, t37 * t41 + t42 * t39 + 0; t40, 0, t38, t40 * pkin(3) + t38 * qJ(4) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:59
	% EndTime: 2020-11-04 20:17:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (23->16), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->8)
	t44 = sin(qJ(3));
	t46 = cos(qJ(3));
	t48 = pkin(3) + pkin(4);
	t49 = t46 * qJ(4) - t48 * t44 - qJ(2);
	t47 = cos(qJ(1));
	t45 = sin(qJ(1));
	t43 = -qJ(5) + pkin(1) + pkin(6);
	t1 = [t45 * t44, -t45 * t46, -t47, t43 * t47 - t49 * t45 + 0; -t47 * t44, t47 * t46, -t45, t43 * t45 + t49 * t47 + 0; t46, t44, 0, t44 * qJ(4) + t48 * t46 + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end