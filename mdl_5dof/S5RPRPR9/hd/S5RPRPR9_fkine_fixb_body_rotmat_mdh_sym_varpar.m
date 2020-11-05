% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:20
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:20:27
	% EndTime: 2020-11-04 20:20:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:20:27
	% EndTime: 2020-11-04 20:20:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [t46, -t45, 0, 0; t45, t46, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:20:27
	% EndTime: 2020-11-04 20:20:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t49 = qJ(1) + pkin(8);
	t48 = cos(t49);
	t47 = sin(t49);
	t1 = [t48, -t47, 0, cos(qJ(1)) * pkin(1) + 0; t47, t48, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:20:27
	% EndTime: 2020-11-04 20:20:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t54 = cos(qJ(3));
	t53 = sin(qJ(3));
	t52 = qJ(1) + pkin(8);
	t51 = cos(t52);
	t50 = sin(t52);
	t1 = [t51 * t54, -t51 * t53, t50, t51 * pkin(2) + t50 * pkin(6) + cos(qJ(1)) * pkin(1) + 0; t50 * t54, -t50 * t53, -t51, t50 * pkin(2) - t51 * pkin(6) + sin(qJ(1)) * pkin(1) + 0; t53, t54, 0, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:20:27
	% EndTime: 2020-11-04 20:20:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->7)
	t58 = sin(qJ(3));
	t59 = cos(qJ(3));
	t60 = pkin(3) * t59 + qJ(4) * t58 + pkin(2);
	t57 = qJ(1) + pkin(8);
	t56 = cos(t57);
	t55 = sin(t57);
	t1 = [t55, -t56 * t59, t56 * t58, cos(qJ(1)) * pkin(1) + t55 * pkin(6) + 0 + t60 * t56; -t56, -t55 * t59, t55 * t58, sin(qJ(1)) * pkin(1) - t56 * pkin(6) + 0 + t60 * t55; 0, -t58, -t59, t58 * pkin(3) - t59 * qJ(4) + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:20:27
	% EndTime: 2020-11-04 20:20:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (47->21), mult. (39->24), div. (0->0), fcn. (52->8), ass. (0->13)
	t72 = pkin(3) + pkin(7);
	t71 = pkin(4) + pkin(6);
	t64 = sin(qJ(5));
	t65 = sin(qJ(3));
	t70 = t64 * t65;
	t66 = cos(qJ(5));
	t69 = t65 * t66;
	t67 = cos(qJ(3));
	t68 = qJ(4) * t65 + t72 * t67 + pkin(2);
	t63 = qJ(1) + pkin(8);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t61 * t66 + t62 * t70, -t61 * t64 + t62 * t69, t62 * t67, cos(qJ(1)) * pkin(1) + 0 + t71 * t61 + t68 * t62; t61 * t70 - t62 * t66, t61 * t69 + t62 * t64, t61 * t67, sin(qJ(1)) * pkin(1) + 0 - t71 * t62 + t68 * t61; -t67 * t64, -t67 * t66, t65, -t67 * qJ(4) + t72 * t65 + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end