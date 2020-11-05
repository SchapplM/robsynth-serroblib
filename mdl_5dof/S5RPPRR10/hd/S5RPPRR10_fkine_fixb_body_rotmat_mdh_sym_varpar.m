% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:16
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:07
	% EndTime: 2020-11-04 20:16:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:07
	% EndTime: 2020-11-04 20:16:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t1 = [t47, -t46, 0, 0; t46, t47, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:07
	% EndTime: 2020-11-04 20:16:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t49 = cos(pkin(8));
	t48 = sin(pkin(8));
	t1 = [t51 * t49, -t51 * t48, t50, t51 * pkin(1) + t50 * qJ(2) + 0; t50 * t49, -t50 * t48, -t51, t50 * pkin(1) - t51 * qJ(2) + 0; t48, t49, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:07
	% EndTime: 2020-11-04 20:16:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t54 = cos(pkin(8));
	t53 = sin(pkin(8));
	t52 = pkin(2) * t54 + qJ(3) * t53 + pkin(1);
	t1 = [t56 * t54, t55, t56 * t53, qJ(2) * t55 + t52 * t56 + 0; t55 * t54, -t56, t55 * t53, -qJ(2) * t56 + t52 * t55 + 0; t53, 0, -t54, pkin(2) * t53 - qJ(3) * t54 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:07
	% EndTime: 2020-11-04 20:16:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t69 = cos(qJ(4));
	t68 = sin(qJ(4));
	t67 = pkin(2) + pkin(3);
	t66 = cos(qJ(1));
	t65 = sin(qJ(1));
	t64 = pkin(6) - qJ(2);
	t63 = cos(pkin(8));
	t62 = sin(pkin(8));
	t59 = t62 * qJ(3) + t67 * t63 + pkin(1);
	t58 = t62 * t69 - t63 * t68;
	t57 = -t62 * t68 - t63 * t69;
	t1 = [-t66 * t57, t66 * t58, -t65, t59 * t66 - t64 * t65 + 0; -t65 * t57, t65 * t58, t66, t59 * t65 + t64 * t66 + 0; t58, t57, 0, -t63 * qJ(3) + t67 * t62 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:07
	% EndTime: 2020-11-04 20:16:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->20), mult. (38->22), div. (0->0), fcn. (52->8), ass. (0->15)
	t76 = qJ(4) + qJ(5);
	t84 = sin(t76);
	t83 = pkin(2) + pkin(3);
	t82 = cos(qJ(1));
	t81 = cos(qJ(4));
	t80 = sin(qJ(1));
	t79 = sin(qJ(4));
	t78 = cos(pkin(8));
	t77 = sin(pkin(8));
	t75 = qJ(2) - pkin(6) - pkin(7);
	t74 = cos(t76);
	t72 = t77 * t74 - t78 * t84;
	t71 = t78 * t74 + t77 * t84;
	t70 = t77 * qJ(3) + t83 * t78 + pkin(1) + (t77 * t79 + t78 * t81) * pkin(4);
	t1 = [t82 * t71, t82 * t72, -t80, t70 * t82 + t75 * t80 + 0; t80 * t71, t80 * t72, t82, t70 * t80 - t75 * t82 + 0; t72, -t71, 0, -t78 * qJ(3) + t83 * t77 + pkin(5) + 0 + (t77 * t81 - t78 * t79) * pkin(4); 0, 0, 0, 1;];
	Tc_mdh = t1;
end