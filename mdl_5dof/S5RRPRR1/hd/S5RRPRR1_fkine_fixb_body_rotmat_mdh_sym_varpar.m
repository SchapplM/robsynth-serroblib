% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:35
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:28
	% EndTime: 2020-11-04 20:35:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:28
	% EndTime: 2020-11-04 20:35:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t1 = [t39, -t38, 0, 0; t38, t39, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:28
	% EndTime: 2020-11-04 20:35:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (4->4), div. (0->0), fcn. (12->4), ass. (0->5)
	t43 = cos(qJ(1));
	t42 = cos(qJ(2));
	t41 = sin(qJ(1));
	t40 = sin(qJ(2));
	t1 = [t43 * t42, -t43 * t40, t41, 0; t41 * t42, -t41 * t40, -t43, 0; t40, t42, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:28
	% EndTime: 2020-11-04 20:35:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (11->9), div. (0->0), fcn. (19->4), ass. (0->7)
	t45 = sin(qJ(1));
	t46 = cos(qJ(2));
	t49 = t45 * t46;
	t47 = cos(qJ(1));
	t48 = t47 * t46;
	t44 = sin(qJ(2));
	t1 = [t48, -t47 * t44, t45, pkin(1) * t48 + t45 * qJ(3) + 0; t49, -t45 * t44, -t47, pkin(1) * t49 - t47 * qJ(3) + 0; t44, t46, 0, t44 * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:28
	% EndTime: 2020-11-04 20:35:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->11), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->9)
	t57 = pkin(2) + pkin(1);
	t58 = cos(qJ(2)) * t57;
	t56 = cos(qJ(1));
	t54 = sin(qJ(1));
	t53 = pkin(3) + qJ(3);
	t52 = qJ(2) + qJ(4);
	t51 = cos(t52);
	t50 = sin(t52);
	t1 = [t56 * t51, -t56 * t50, t54, t53 * t54 + t56 * t58 + 0; t54 * t51, -t54 * t50, -t56, -t56 * t53 + t54 * t58 + 0; t50, t51, 0, sin(qJ(2)) * t57 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:28
	% EndTime: 2020-11-04 20:35:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->16), mult. (28->20), div. (0->0), fcn. (41->8), ass. (0->15)
	t63 = sin(qJ(5));
	t64 = sin(qJ(1));
	t73 = t64 * t63;
	t65 = cos(qJ(5));
	t72 = t64 * t65;
	t67 = cos(qJ(1));
	t71 = t67 * t63;
	t70 = t67 * t65;
	t61 = qJ(2) + qJ(4);
	t59 = sin(t61);
	t68 = pkin(2) + pkin(1);
	t69 = pkin(4) * t59 + cos(qJ(2)) * t68;
	t62 = pkin(3) + qJ(3);
	t60 = cos(t61);
	t1 = [t60 * t70 + t73, -t60 * t71 + t72, t67 * t59, t62 * t64 + t69 * t67 + 0; t60 * t72 - t71, -t60 * t73 - t70, t64 * t59, -t67 * t62 + t69 * t64 + 0; t59 * t65, -t59 * t63, -t60, -t60 * pkin(4) + sin(qJ(2)) * t68 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end