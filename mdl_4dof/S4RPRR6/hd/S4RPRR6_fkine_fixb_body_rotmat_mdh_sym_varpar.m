% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:44
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:15
	% EndTime: 2020-11-04 19:44:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:15
	% EndTime: 2020-11-04 19:44:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t28 = cos(qJ(1));
	t27 = sin(qJ(1));
	t1 = [t28, -t27, 0, 0; t27, t28, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:15
	% EndTime: 2020-11-04 19:44:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t32 = cos(qJ(1));
	t31 = sin(qJ(1));
	t30 = cos(pkin(7));
	t29 = sin(pkin(7));
	t1 = [t32 * t30, -t32 * t29, t31, t32 * pkin(1) + t31 * qJ(2) + 0; t31 * t30, -t31 * t29, -t32, t31 * pkin(1) - t32 * qJ(2) + 0; t29, t30, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:15
	% EndTime: 2020-11-04 19:44:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t37 = pkin(5) + qJ(2);
	t36 = pkin(7) + qJ(3);
	t35 = cos(t36);
	t34 = sin(t36);
	t33 = cos(pkin(7)) * pkin(2) + pkin(1);
	t1 = [t39 * t35, -t39 * t34, t38, t39 * t33 + t37 * t38 + 0; t38 * t35, -t38 * t34, -t39, t38 * t33 - t39 * t37 + 0; t34, t35, 0, sin(pkin(7)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:15
	% EndTime: 2020-11-04 19:44:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t45 = pkin(7) + qJ(3);
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t44 = -pkin(6) - pkin(5) - qJ(2);
	t43 = qJ(4) + t45;
	t42 = cos(t43);
	t41 = sin(t43);
	t40 = pkin(3) * cos(t45) + cos(pkin(7)) * pkin(2) + pkin(1);
	t1 = [t47 * t42, -t47 * t41, t46, t40 * t47 - t44 * t46 + 0; t46 * t42, -t46 * t41, -t47, t40 * t46 + t44 * t47 + 0; t41, t42, 0, pkin(3) * sin(t45) + sin(pkin(7)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end