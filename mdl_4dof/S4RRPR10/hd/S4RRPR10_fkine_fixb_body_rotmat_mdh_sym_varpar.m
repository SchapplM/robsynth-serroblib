% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:48:28
	% EndTime: 2020-11-04 19:48:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:48:28
	% EndTime: 2020-11-04 19:48:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = cos(qJ(1));
	t30 = sin(qJ(1));
	t1 = [t31, -t30, 0, 0; t30, t31, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:48:28
	% EndTime: 2020-11-04 19:48:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t35 = cos(qJ(1));
	t34 = cos(qJ(2));
	t33 = sin(qJ(1));
	t32 = sin(qJ(2));
	t1 = [t35 * t34, -t35 * t32, t33, t35 * pkin(1) + t33 * pkin(5) + 0; t33 * t34, -t33 * t32, -t35, t33 * pkin(1) - t35 * pkin(5) + 0; t32, t34, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:48:28
	% EndTime: 2020-11-04 19:48:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t40 = cos(qJ(1));
	t39 = cos(qJ(2));
	t38 = sin(qJ(1));
	t37 = sin(qJ(2));
	t36 = t39 * pkin(2) + t37 * qJ(3) + pkin(1);
	t1 = [t38, -t40 * t39, t40 * t37, t38 * pkin(5) + t36 * t40 + 0; -t40, -t38 * t39, t38 * t37, -t40 * pkin(5) + t36 * t38 + 0; 0, -t37, -t39, t37 * pkin(2) - t39 * qJ(3) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:48:28
	% EndTime: 2020-11-04 19:48:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t42 = sin(qJ(4));
	t44 = sin(qJ(1));
	t53 = t44 * t42;
	t45 = cos(qJ(4));
	t52 = t44 * t45;
	t47 = cos(qJ(1));
	t51 = t47 * t42;
	t50 = t47 * t45;
	t49 = pkin(2) + pkin(6);
	t48 = pkin(3) + pkin(5);
	t46 = cos(qJ(2));
	t43 = sin(qJ(2));
	t41 = qJ(3) * t43 + t46 * t49 + pkin(1);
	t1 = [t43 * t51 + t52, t43 * t50 - t53, t47 * t46, t41 * t47 + t44 * t48 + 0; t43 * t53 - t50, t43 * t52 + t51, t44 * t46, t41 * t44 - t47 * t48 + 0; -t46 * t42, -t46 * t45, t43, -qJ(3) * t46 + t43 * t49 + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end