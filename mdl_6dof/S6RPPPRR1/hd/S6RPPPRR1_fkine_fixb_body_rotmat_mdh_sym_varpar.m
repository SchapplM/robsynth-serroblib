% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPPRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:22
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:53
	% EndTime: 2020-11-04 21:22:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:53
	% EndTime: 2020-11-04 21:22:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t43 = cos(qJ(1));
	t42 = sin(qJ(1));
	t1 = [t43, -t42, 0, 0; t42, t43, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:53
	% EndTime: 2020-11-04 21:22:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t46 = qJ(1) + pkin(9);
	t45 = cos(t46);
	t44 = sin(t46);
	t1 = [t45, -t44, 0, cos(qJ(1)) * pkin(1) + 0; t44, t45, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:53
	% EndTime: 2020-11-04 21:22:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t49 = qJ(1) + pkin(9);
	t48 = cos(t49);
	t47 = sin(t49);
	t1 = [0, -t48, t47, t48 * pkin(2) + t47 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; 0, -t47, -t48, t47 * pkin(2) - t48 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; 1, 0, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:53
	% EndTime: 2020-11-04 21:22:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (20->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t53 = pkin(2) + qJ(4);
	t52 = qJ(1) + pkin(9);
	t51 = cos(t52);
	t50 = sin(t52);
	t1 = [0, t50, t51, t53 * t51 + t50 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; 0, -t51, t50, t53 * t50 - t51 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; 1, 0, 0, pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:53
	% EndTime: 2020-11-04 21:22:53
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->15), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->8)
	t60 = cos(qJ(5));
	t59 = sin(qJ(5));
	t58 = pkin(2) + qJ(4);
	t57 = pkin(7) - qJ(3);
	t56 = qJ(1) + pkin(9);
	t55 = cos(t56);
	t54 = sin(t56);
	t1 = [t55 * t59, t55 * t60, -t54, t58 * t55 - t57 * t54 + cos(qJ(1)) * pkin(1) + 0; t54 * t59, t54 * t60, t55, t57 * t55 + t58 * t54 + sin(qJ(1)) * pkin(1) + 0; t60, -t59, 0, pkin(4) + pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:53
	% EndTime: 2020-11-04 21:22:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->24), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->12)
	t66 = sin(qJ(6));
	t67 = sin(qJ(5));
	t72 = t66 * t67;
	t68 = cos(qJ(6));
	t71 = t67 * t68;
	t69 = cos(qJ(5));
	t70 = pkin(5) * t67 - pkin(8) * t69 + pkin(2) + qJ(4);
	t64 = pkin(7) - qJ(3);
	t63 = qJ(1) + pkin(9);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [-t61 * t66 + t62 * t71, -t61 * t68 - t62 * t72, -t62 * t69, cos(qJ(1)) * pkin(1) - t64 * t61 + 0 + t70 * t62; t61 * t71 + t62 * t66, -t61 * t72 + t62 * t68, -t61 * t69, sin(qJ(1)) * pkin(1) + t64 * t62 + 0 + t70 * t61; t69 * t68, -t69 * t66, t67, t69 * pkin(5) + t67 * pkin(8) + pkin(3) + pkin(4) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end