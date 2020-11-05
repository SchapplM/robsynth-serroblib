% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR14 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:21
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:49
	% EndTime: 2020-11-04 20:21:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:49
	% EndTime: 2020-11-04 20:21:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [t44, -t43, 0, 0; t43, t44, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:49
	% EndTime: 2020-11-04 20:21:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [0, -t46, t45, t46 * pkin(1) + t45 * qJ(2) + 0; 0, -t45, -t46, t45 * pkin(1) - t46 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:49
	% EndTime: 2020-11-04 20:21:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t51 = pkin(1) + pkin(6);
	t50 = cos(qJ(1));
	t49 = cos(qJ(3));
	t48 = sin(qJ(1));
	t47 = sin(qJ(3));
	t1 = [t48 * t47, t48 * t49, t50, t48 * qJ(2) + t51 * t50 + 0; -t50 * t47, -t50 * t49, t48, -t50 * qJ(2) + t51 * t48 + 0; t49, -t47, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:49
	% EndTime: 2020-11-04 20:21:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t56 = qJ(3) + pkin(8);
	t55 = pkin(1) + pkin(6) + qJ(4);
	t54 = cos(t56);
	t53 = sin(t56);
	t52 = sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t57 * t53, t57 * t54, t58, t52 * t57 + t55 * t58 + 0; -t58 * t53, -t58 * t54, t57, -t52 * t58 + t55 * t57 + 0; t54, -t53, 0, cos(qJ(3)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:49
	% EndTime: 2020-11-04 20:21:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->22), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t67 = sin(qJ(5));
	t69 = sin(qJ(1));
	t77 = t69 * t67;
	t70 = cos(qJ(5));
	t76 = t69 * t70;
	t72 = cos(qJ(1));
	t75 = t72 * t67;
	t74 = t72 * t70;
	t65 = sin(pkin(8));
	t66 = cos(pkin(8));
	t71 = cos(qJ(3));
	t73 = (t66 * pkin(4) + t65 * pkin(7) + pkin(3)) * sin(qJ(3)) - (-t65 * pkin(4) + t66 * pkin(7)) * t71 + qJ(2);
	t64 = qJ(3) + pkin(8);
	t63 = pkin(1) + pkin(6) + qJ(4);
	t62 = cos(t64);
	t61 = sin(t64);
	t1 = [t61 * t76 + t75, -t61 * t77 + t74, -t69 * t62, t63 * t72 + t73 * t69 + 0; -t61 * t74 + t77, t61 * t75 + t76, t72 * t62, t63 * t69 - t73 * t72 + 0; t62 * t70, -t62 * t67, t61, t71 * pkin(3) + t62 * pkin(4) + t61 * pkin(7) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end