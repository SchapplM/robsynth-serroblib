% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:32
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:34
	% EndTime: 2020-11-04 19:32:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:34
	% EndTime: 2020-11-04 19:32:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = cos(pkin(6));
	t30 = sin(pkin(6));
	t1 = [t31, -t30, 0, 0; t30, t31, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:34
	% EndTime: 2020-11-04 19:32:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t33 = cos(pkin(6));
	t32 = sin(pkin(6));
	t1 = [t33, 0, t32, t33 * pkin(1) + t32 * qJ(2) + 0; t32, 0, -t33, t32 * pkin(1) - t33 * qJ(2) + 0; 0, 1, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:34
	% EndTime: 2020-11-04 19:32:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t40 = pkin(1) + pkin(2);
	t39 = cos(qJ(3));
	t38 = sin(qJ(3));
	t37 = cos(pkin(6));
	t36 = sin(pkin(6));
	t35 = t36 * t39 - t37 * t38;
	t34 = -t36 * t38 - t37 * t39;
	t1 = [-t34, t35, 0, t36 * qJ(2) + t40 * t37 + 0; t35, t34, 0, -t37 * qJ(2) + t40 * t36 + 0; 0, 0, -1, -pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:34
	% EndTime: 2020-11-04 19:32:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->19), mult. (32->20), div. (0->0), fcn. (46->6), ass. (0->12)
	t51 = pkin(1) + pkin(2);
	t50 = cos(qJ(3));
	t49 = cos(qJ(4));
	t48 = sin(qJ(3));
	t47 = sin(qJ(4));
	t46 = cos(pkin(6));
	t45 = sin(pkin(6));
	t44 = t46 * pkin(3) - pkin(5) * t45;
	t43 = pkin(3) * t45 + t46 * pkin(5);
	t42 = t45 * t48 + t46 * t50;
	t41 = -t45 * t50 + t46 * t48;
	t1 = [t42 * t49, -t42 * t47, t41, t45 * qJ(2) + t43 * t48 + t44 * t50 + t51 * t46 + 0; -t41 * t49, t41 * t47, t42, -t46 * qJ(2) + t43 * t50 - t44 * t48 + t51 * t45 + 0; -t47, -t49, 0, -pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end